#ifndef FPS_TOPOLOGY_LOCAL_GAUSS
#define FPS_TOPOLOGY_LOCAL_GAUSS

#include "fstream"
#include <sstream>        // std::stringstream, std::stringbuf
#include <iostream>       // std::cout, std::fixed
#include <iomanip>        // std::setprecision
#include <ctime>          // std::time
#include <algorithm>      // std::remove_if
#include "stdio.h"
#include <cmath>
#include <vector>
#include <list>
#include <random>
#include "assert.h"
#include <memory>         // std::unique_ptr

#include "rng.h"             // random numbers, distributions
#include "helper_io.hpp"     // hdf5
#include "helper_types.hpp"  // neuron, electrode, topology

// choose the precision for custom float, relevant for stuff to write to disk
typedef float real_t;
#define H5T_CUSTOM_REAL (H5T_NATIVE_FLOAT)
// typedef double real_t;
// #define H5T_CUSTOM_REAL (H5T_NATIVE_DOUBLE)

// save some space by aliasing the vector class
template<class T1>
using vec = std::vector <T1>;

class topology_local_gauss {
 public:

  // global parameter, we need useful setters that alter values consistently
  struct parameters {
    size_t N   = 160000;                       // number neurons
    size_t K   = 1000;                         // number outgoing connection
    double rho = 100.0;                        // [neurons/mm2]
    double L   = 1000.0*sqrt(double(N)/rho);   // [um] dish size
    double std = 6.0;                          // [d_N] eff. conn-length

    // choose d so that approx K=num_outgoing neighbours are accessible
    double d_max = sqrt(L *L *K/M_PI/double(N)); // [um] max dist. for connect

    // analytic solution (ideal gas) for nearest-neuron distance d_N
    double d_N = 0.5*L/sqrt(double(N));        // [um]
  } par;

  // adjusts L, rho, d_max accordingly
  void set_N_and_d_N(size_t N, double d_N) {
    par.N     = N;
    par.d_N   = d_N;
    par.L     = 2.0*sqrt(double(par.N))*par.d_N;
    par.rho   = double(par.N)/pow(par.L/1000.0, 2.0);
    par.d_max = sqrt(par.L*par.L*par.K/M_PI/double(par.N));
  }

  void print_parameters() {
    printf("local gauss topology prameters:\n");
    printf("\tN       %lu [neurons]\n", par.N);
    printf("\tK       %lu [neurons]\n", par.K);
    printf("\trho     %.2e [neurons/mm2]\n", par.rho);
    printf("\tL       %.2e [um]\n", par.L);
    printf("\tstd     %.2e [d_N] = %.2e [um]\n", par.std, par.std*par.d_N);
    printf("\td_max   %.2e [um]\n", par.d_max);
    printf("\td_N     %.2e [um]\n", par.d_N);
  }

  // find distance between neurons (or electrodes)
  template <typename T1, typename T2>
  inline double get_dist_squ(T1 *p, T2 *t) {
    double xr = std::max(t->x,p->x);
    double xl = std::min(t->x,p->x);
    double yr = std::max(t->y,p->y);
    double yl = std::min(t->y,p->y);

    // pbc
    double x2 = std::min((xr-xl)*(xr-xl), (par.L-xr+xl)*(par.L-xr+xl));
    double y2 = std::min((yr-yl)*(yr-yl), (par.L-yr+yl)*(par.L-yr+yl));

    return (x2+y2);
  }

  // ------------------------------------------------------------------ //
  // create topology from parameters set above.
  // hid_t is just an int.
  // deadzones where not to place neurons can be passed using the
  // avoid_ vector.
  // ------------------------------------------------------------------ //

  template <typename T1>
  void init(vec <neuron *> &neurons, hid_t h5file, const vec <T1 *> &avoid_xy, double avoid_radius_squ) {

    // need to add a check that it is possible to distribute all the neurons
    // in the given space
    printf("placing neurons on local gauss topology\n");
    for (size_t i = 0; i < neurons.size(); i++) delete neurons[i];
    neurons.resize(0);
    neurons.reserve(par.N);
    set_ndis_param(0.0, 1.0);   // normal dist mu=0 std=1
    set_udis_param(0.0, 1.0);   // uniform dist min=0 max=1

    // create and place neurons
    size_t rejections = 0;
    for (size_t i = 0; i < par.N; i++) {
      neuron *n = new neuron(i);
      bool placed = false;
      while (!placed) {
        n->x = udis(rng)*par.L;
        n->y = udis(rng)*par.L;
        placed = true;
        for (size_t a = 0; a < avoid_xy.size(); a++) {
          if (get_dist_squ(avoid_xy[a], n) <= avoid_radius_squ) {
            placed = false;
            rejections += 1;
            if (rejections > par.N/10) {
              printf("rejected too many neuron placements (%lu) ", rejections);
              printf("are the parameters set correctly?\n");
              exit(-1);
            }
            break;
          }
        }
      }
      neurons.push_back(n);
    }

    // create physical connections with hard limit given by par.d_max
    printf("connecting neurons with radius %.2e [um]\n", par.d_max);
    double dist_squ_limit  = par.d_max*par.d_max;
    size_t num_connections = 0;
    for (size_t i = 0; i < neurons.size(); i++) {
      neuron *src = neurons[i];
      vec <std::pair <double, neuron *> > sortable;
      sortable.reserve(size_t(par.K*1.1));

      for (size_t j = 0; j < neurons.size(); j++) {
        neuron *tar = neurons[j];
        if (src == tar) continue;
        double dist_squ = get_dist_squ(src, tar);
        if (dist_squ < dist_squ_limit) {
          sortable.push_back( std::make_pair(dist_squ, tar) );
          num_connections += 1;
        }
      }
      // sort outgoing connections by distance
      sort(sortable.begin(), sortable.end());
      src->outgoing.reserve(sortable.size());
      for (size_t j = 0; j < sortable.size(); j++)
        src->outgoing.push_back(sortable[j].second);

      if(i==0 || is_percent(i, size_t(neurons.size()), 10.)) {
        printf("\t%s, %lu/%lu\n", time_now().c_str(), i, neurons.size());
      }
    }

    // set synapctic weight (effective interaction radius) via std (sigma)
    double sigma_squ = 2.*pow(par.std*par.d_N, 2.);
    for (size_t i = 0; i < neurons.size(); i++) {
      neuron *src   = neurons[i];
      size_t n_cout = src->outgoing.size();
      src->outgoing_probability.reserve(n_cout);
      double norm = 0.;

      // activation probabilities based on distance
      for (size_t j = 0; j < n_cout; j++) {
        double dij_squ = get_dist_squ(src, src->outgoing[j]);
        double pij     = exp(-dij_squ/sigma_squ);
        src->outgoing_probability.push_back(pij);
        norm += pij;
      }

      // normalize probabilities
      for (size_t j = 0; j < n_cout; j++) {
        src->outgoing_probability[j] /= norm;
      }
    }

    // print details
    printf("init done\n");
    printf("number of neurons placed: %lu\n", neurons.size());
    printf("\trejections in dead zones: %lu\n", rejections);
    printf("\toutgoing connections per neuron: K=%lu\n", par.K);
    printf("\t\t(measured: %.2f)\n", num_connections/double(neurons.size()));
    printf("\tconnections created within: d_max=%.2e [um]\n", par.d_max);
    printf("\taverage distance between nearest neighbours: %.2e [um]\n", par.d_N);
    printf("\t\t(measured: %e)\n",
           measure_avg_nearest_neighbour_distance(neurons));
  }

  // constructor overload. call init after setting custom parameters
  void init(vec <neuron *> &neurons, hid_t h5file = -1) {
    vec <neuron *> placeholder;
    init(neurons, h5file, placeholder, 0.0);
  }

  // ------------------------------------------------------------------ //
  // debug helper
  // ------------------------------------------------------------------ //

  double measure_avg_nearest_neighbour_distance(vec <neuron *> &neurons) {
    double avg_dist = 0;
    for (size_t i = 0; i < neurons.size(); i++) {
      double shortest_distance = std::numeric_limits<double>::max();
      for (size_t j = 0; j < neurons[i]->outgoing.size(); j++) {
        // if (neurons[i] == j) continue;
        assert(neurons[i] != neurons[i]->outgoing[j]);
        double d = get_dist_squ(neurons[i], neurons[i]->outgoing[j]);
        if (d < shortest_distance) shortest_distance = d;
      }
      avg_dist += sqrt(shortest_distance);
    }
    avg_dist = avg_dist / double(neurons.size());
    return avg_dist;
  }

  double measure_avg_connection_length(vec <neuron *> &neurons) {
    double avg_dist = 0;
    for (size_t i = 0; i < neurons.size(); i++) {
      for (size_t j = 0; j < neurons[i]->outgoing.size(); j++) {
        double d = get_dist_squ(neurons[i], neurons[i]->outgoing[j]);
        avg_dist += sqrt(d);
      }
    }
    avg_dist = avg_dist / double(par.N) / double(par.K);
    return avg_dist;
  }

  // d_N (former name delta_l)
  double measure_avg_distance(vec <neuron *> &neurons) {
    double dl = 0;
    for (size_t i = 0; i < neurons.size(); i++) {
      for (size_t j = 0; j < neurons.size(); j++) {
        double d = get_dist_squ(neurons[i], neurons[j]);
        dl += sqrt(d);
      }
    }
    dl = dl / double(neurons.size()) / double(neurons.size());
    return dl;
  }

  // ------------------------------------------------------------------ //
  // helpers to write details to hd5f
  // h5file should be a file (H5Fcreate) with groups (H5Gcreate):
  // /neurons /meta
  // ------------------------------------------------------------------ //

  void write_topology_details(hid_t h5file) {
    if (h5file < 0) {
      printf("invalid h5file, no topology details written");
      return;
    } else {
      hdf5_write_string(h5file, "/meta/topology", "local gauss");
      hdf5_write_scalar(h5file, "/meta/num_neur",
                        par.N, H5T_NATIVE_HSIZE);
      hdf5_write_scalar(h5file, "/meta/num_outgoing",
                        par.K, H5T_NATIVE_HSIZE);
      hdf5_write_scalar(h5file, "/meta/sys_size",
                        par.L, H5T_NATIVE_DOUBLE);
      hdf5_write_scalar(h5file, "/meta/neuron_density",
                        par.rho, H5T_NATIVE_DOUBLE);
      hdf5_write_scalar(h5file, "/meta/d_max",
                        par.d_max, H5T_NATIVE_DOUBLE);
      hdf5_write_scalar(h5file, "/meta/sigma",
                        par.std, H5T_NATIVE_DOUBLE);
      hdf5_write_scalar(h5file, "/meta/neur_dist",
                        par.d_N, H5T_NATIVE_DOUBLE);
    }
  }

  void write_connectivty_matrix(vec <neuron *> &neurons, hid_t h5file) {
    hid_t h5_matrix = hdf5_create_appendable_nd(
                          h5file, "/connectivity_matrix", H5T_NATIVE_HSIZE,
                          neurons.size(), neurons.size());
    for (size_t i = 0; i < neurons.size(); i++) {
      vec <size_t> row(neurons.size(), 0);
      for (size_t j = 0; j < neurons[i]->outgoing.size(); j++)
        row[ neurons[i]->outgoing[j]->id ] = 1;
      hdf5_append_nd(h5_matrix, row, H5T_NATIVE_HSIZE, i, 0);
    }
    H5Dclose(h5_matrix);
  }

  void write_neuron_details(vec <neuron *> &neurons, hid_t h5file) {
    vec <real_t> x(neurons.size(), 0);
    vec <real_t> y(neurons.size(), 0);
    for (size_t i = 0; i < neurons.size(); i++) {
      x[i] = neurons[i]->x;
      y[i] = neurons[i]->y;
    }

    hdf5_write_vector(h5file, "/neurons/pos_x",
                      x, H5T_CUSTOM_REAL);
    hdf5_write_vector(h5file, "/neurons/pos_y",
                      y, H5T_CUSTOM_REAL);
  }
};

#endif
