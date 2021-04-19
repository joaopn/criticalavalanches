#ifndef FPS_TOPOLOGY_ORLANDI
#define FPS_TOPOLOGY_ORLANDI

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

class axon {
 public:
  double x;
  double y;
  double l;
  double phi;
  size_t ns;
  vec <real_t> x_path;
  vec <real_t> y_path;
};

class topology_orlandi {
 public:

  // global parameter
  struct parameters {
    size_t N   = 1600;                         // number neurons
    double rho = 100.0;                        // [neurons/mm2]
    double L   = 1000.0*sqrt(double(N)/rho);   // [um] dish size

    // axons, variable length, segments of fixed length with variable angle
    double std_l   = 800.0;  // [um]  st. dev. of rayleigh dist. for axon len
    double max_l   = 1500.0; // [um]  max lentgh allowed for axons
    double del_l   = 10.0;   // [um]  length of each segment of axons
    double std_phi = 0.1;    // [rad] std of Gauss dist. for angles betw segs

    // soma, hard sphere
    double R_s     = 7.5;    // [um] radius of soma

    // dendritic tree, sphere with variable diameter
    double std_d   = 20.0;   // [um] std of Gauss dist. for diameter
    double mu_d    = 150.0;  // [um] mean of Gauss dist. for diameter

    // connection probability alpha when axons intersect dendritic tree
    double alpha = .66;

    // only needed for export consistent with other topologies
    double K   = 0;           // set later, number realized outgoing conn.
    double d_N = 1000./2./sqrt(rho);
  } par;

  void print_parameters() {
    printf("orlandi topology parameters:\n");
    printf("\tN       %lu  [neurons] \n", par.N);
    printf("\trho     %.2e [neurons/mm2]\n", par.rho);
    printf("\tL       %.2e [um]\n", par.L);
    printf("\tstd_l   %.2e [um]\n", par.std_l);
    printf("\tmax_l   %.2e [um]\n", par.max_l);
    printf("\tdel_l   %.2e [um]\n", par.del_l);
    printf("\tstd_phi %.2e [rad]\n", par.std_phi);
    printf("\tR_s     %.2e [um]\n", par.R_s);
    printf("\tstd_d   %.2e [um]\n", par.std_d);
    printf("\tmu_d    %.2e [um]\n", par.mu_d);
    printf("\talpha   %.2e\n", par.alpha);
  }

  // find distance between neurons (or electrodes)
  template <typename T1, typename T2>
  inline double get_dist_squ(T1 *p, T2 *t) {
    double xr = std::max(t->x,p->x);
    double xl = std::min(t->x,p->x);
    double yr = std::max(t->y,p->y);
    double yl = std::min(t->y,p->y);

    // pbc
    // double x2 = std::min((xr-xl)*(xr-xl), (par.L-xr+xl)*(par.L-xr+xl));
    // double y2 = std::min((yr-yl)*(yr-yl), (par.L-yr+yl)*(par.L-yr+yl));

    double x2 = (xr-xl)*(xr-xl);
    double y2 = (yr-yl)*(yr-yl);

    return (x2+y2);
  }

  // ------------------------------------------------------------------ //
  // domain decomposition. we need to compare every axon segment
  // with every dendritic tree. this should make it a bit faster.
  // ------------------------------------------------------------------ //

  double box_size;
  int nb;
  vec <vec <neuron *> > boxes;
  vec <vec <int> >      boxes_to_check;

  inline size_t box_id_from_pos(double x, double y) {
    int ix     = std::floor(x/box_size)+1;
    int iy     = std::floor(y/box_size)+1;
    size_t bid = size_t(ix + nb*iy);

    assert(x   < par.L+box_size);
    assert(y   < par.L+box_size);
    assert(x   > 0 -box_size);
    assert(y   > 0 -box_size);
    assert(bid < boxes.size());

    return bid;
  }

  void init_domain_decomposition(double min_box_size) {
    printf("domain decomposition\n");

    // nb boxes per dim, 2x border padding so we dont cut acons immediately
    nb       = size_t((par.L)/min_box_size)+3;
    box_size = double(par.L + 2.0*min_box_size)/double(nb);
    boxes_to_check.resize(nb*nb);
    boxes.resize(nb*nb);

    for (int ix = 0; ix < nb; ix++) {
      for (int iy = 0; iy < nb; iy++) {
        int i = ix + nb*iy;
        vec <int> tar;
        tar = {i, i-1, i+1, i-nb, i+nb, i-nb-1, i-nb+1, i+nb-1, i+nb+1};

        // remove guys that are not wanted (boundaries)
        vec <int> ex;
        if      (ix == 0)    ex.insert(ex.end(), {i-1,  i-nb-1, i+nb-1});
        else if (ix == nb-1) ex.insert(ex.end(), {i+1,  i-nb+1, i+nb+1});
        if      (iy == 0)    ex.insert(ex.end(), {i-nb, i-nb-1, i-nb+1});
        else if (iy == nb-1) ex.insert(ex.end(), {i+nb, i+nb-1, i+nb+1});

        for (size_t j = 0; j < ex.size(); j++)
          tar.erase(std::remove(tar.begin(), tar.end(), ex[j]), tar.end());

        boxes_to_check[i] = tar;
        boxes[i] = vec <neuron *> ();
      }
    }
  }

  // ------------------------------------------------------------------ //
  // create topology from parameters set above.
  // hid_t is just an int. if < 0 nothing will be written, else the axon
  // topology will be written to the specified hdf5 file.
  // deadzones where not to place neurons can be passed using the
  // avoid_ vector
  // ------------------------------------------------------------------ //

  template <typename T1>
  void init(vec <neuron *> &neurons, hid_t h5file, const vec <T1 *> &avoid_xy, double avoid_radius_squ) {

    // need to add a check that it is possible to distribute all the neurons
    // in the given space
    printf("placing neurons on orlandi topology\n");
    size_t num_outgoing = 0;
    size_t num_neur     = 0;

    for (size_t i = 0; i < neurons.size(); i++) delete neurons[i];
    neurons.resize(0);
    neurons.reserve(par.N);

    // gaussian distribution for dendritic tree size
    set_ndis_param(par.mu_d, par.std_d);
    double max_tree_radius = 0;
    double min_dist_squ    = pow(2.0*par.R_s, 2.0);

    // draw dendritic radii, create neurons with position outside the domain
    for (size_t i = 0; i < par.N; i++) {
      neuron *n = new neuron(i);
      n->R_d = ndis(rng);
      if (n->R_d > max_tree_radius) max_tree_radius = n->R_d;
      neurons.push_back(n);
    }

    init_domain_decomposition(max_tree_radius);

    // place neurons in domain decomposition boxes
    for (size_t i = 0; i < neurons.size(); i++) {
      neuron *n = neurons[i];
      bool placed = false;
      size_t bid;
      while (!placed) {
        placed = true;
        n->x = udis(rng)*par.L;
        n->y = udis(rng)*par.L;
        for (size_t a = 0; a < avoid_xy.size(); a++) {
          if (get_dist_squ(avoid_xy[a], n) <= avoid_radius_squ)
            placed = false;
        }
        bid  = box_id_from_pos(n->x, n->y);
        for (size_t b = 0; b < boxes_to_check[bid].size(); b++) {
          size_t btc = boxes_to_check[bid][b];
          for (size_t bb = 0; bb < boxes[btc].size(); bb++) {
            size_t j = boxes[btc][bb]->id;
            if (get_dist_squ(neurons[j], n) <= min_dist_squ)
              placed = false;
          }
        }
      }
      boxes[bid].push_back(n);
    }

    // print an overview of the domain decomposition
    printf("neurons per box (%d * %d boxes)\n", nb, nb);
    for (int ix = 0; ix < nb; ix++) {
      for (int iy = 0; iy < nb; iy++) {
        int i = ix + nb*iy;
        num_neur += boxes[i].size();
        // printf("| %*d: %*lu ", 3, i, 3, boxes[i].size());
        // printf("| %*lu ", 2, boxes[i].size());
      }
      // printf("|\n");
    }
    assert(num_neur == neurons.size());

    // ----------------------------------------------------------------- //
    // grow axons and create connections
    // ----------------------------------------------------------------- //

    hid_t h5_x, h5_y;
    if (h5file >= 0) {
      // create hdf5 data sets
      h5_x = hdf5_create_appendable_nd(
                 h5file, "/axons/segments_x", H5T_CUSTOM_REAL,
                 neurons.size(), 1000);
      h5_y = hdf5_create_appendable_nd(
                 h5file, "/axons/segments_y", H5T_CUSTOM_REAL,
                 neurons.size(), 1000);
    }

    set_ndis_param(0.0, 1.0);   // normal dist mu=0 std=1
    set_udis_param(0.0, 1.0);   // uniform dist min=0 max=1
    set_rdis_param(par.std_l);  // rayleigh dist with std (=800)

    for (size_t i = 0; i < neurons.size(); i++) {
      neuron *n = neurons[i];
      axon  *ax = new axon();
      ax->l     = rdis(rng);                          // total axon length
      ax->phi   = udis(rng)*2.0*M_PI;                 // [rad] init direction
      ax->x     = n->x + par.R_s*cos(ax->phi);        // [um]  init position
      ax->y     = n->y + par.R_s*sin(ax->phi);
      ax->ns    = size_t(ax->l/par.del_l)+1;          // number segments
      n ->l     = ax->l;

      // avoid N*N memory requirement and provide fast check for crossing.
      vec <bool> intersec(neurons.size(), false);
      size_t num_intersec = 0;

      // construct chain of segments
      for (size_t s = 0; s < ax->ns; s++) {
        if (h5file >= 0) {
          ax->x_path.push_back(ax->x);
          ax->y_path.push_back(ax->y);
        }
        ax->phi += par.std_phi*ndis(rng);
        ax->x   += par.del_l*cos(ax->phi);
        ax->y   += par.del_l*sin(ax->phi);

        // stop growing when axons grows out of the substrate
        if (ax->x >= par.L+box_size || ax->y >= par.L+box_size ||
            ax->x <= 0    -box_size || ax->y <= 0    -box_size)
          break;

        // find box of axon segment
        size_t bid = box_id_from_pos(ax->x, ax->y);
        for (size_t b = 0; b < boxes_to_check[bid].size(); b++) {
          size_t btc = boxes_to_check[bid][b];
          for (size_t bb = 0; bb < boxes[btc].size(); bb++) {
            size_t j = boxes[btc][bb]->id;
            if (!intersec[j]) {
              neuron *m = neurons[j];
              if (get_dist_squ(ax, m) <= m->R_d*m->R_d) {
                intersec[j] = true;
                num_intersec += 1;
              }
            }
          }
        }
      }

      // write axons segment details if output file available
      if (h5file >= 0) {
        ax->x_path.push_back(ax->x);
        ax->y_path.push_back(ax->y);
        hdf5_append_nd(h5_x, ax->x_path, H5T_CUSTOM_REAL, i, 0);
        hdf5_append_nd(h5_y, ax->y_path, H5T_CUSTOM_REAL, i, 0);
      }

      // write pointer list for outgoing connections. corresponds to
      // connection matrix.
      // sort list by distance between neurons for coalesence compensation
      vec <std::pair <double, neuron *> > sortable;
      for (size_t j = 0; j < intersec.size(); j++) {
        neuron *tar = neurons[j];
        if (intersec[j] && n != tar) {
          if (udis(rng) < par.alpha) {
            double dist_squ = get_dist_squ(n, tar);
            sortable.push_back( std::make_pair(dist_squ, tar) );
            num_outgoing += 1;
          } else
            // remember failed connection attempts  so we can plot them.
            n->rejected.push_back(neurons[j]);
        }
      }

      // sort outgoing connections by distance
      sort(sortable.begin(), sortable.end());
      for (size_t j = 0; j < sortable.size(); j++)
        n->outgoing.push_back(sortable[j].second);

      // set (constant) outgoing porbablity to 1/number_connections
      if (sortable.size() > 0)
        n->outgoing_probability.push_back(1./double(sortable.size()));

      delete ax;
    }

    par.K = double(num_outgoing)/double(neurons.size());

    if (h5file >= 0) {
      H5Dclose(h5_x);
      H5Dclose(h5_y);
    }
  }

  // constructor overload. call init after setting custom parameters
  void init(vec <neuron *> &neurons, hid_t h5file = -1) {
    vec <neuron *> placeholder;
    init(neurons, h5file, placeholder, 0.0);
  }

  // ------------------------------------------------------------------ //
  // helpers to write details to hd5f
  // h5file should be a file (H5Fcreate) with groups (H5Gcreate):
  // /neurons and /axons
  // ------------------------------------------------------------------ //

  void write_topology_details(hid_t h5file) {
    if (h5file < 0) {
      printf("invalid h5file, no topology details written");
      return;
    } else {
      hdf5_write_string(h5file, "/meta/topology", "orlandi");
      hdf5_write_scalar(h5file, "/meta/num_neur",
                        par.N, H5T_NATIVE_HSIZE);
      hdf5_write_scalar(h5file, "/meta/num_outgoing",
                        par.K, H5T_NATIVE_DOUBLE);
      hdf5_write_scalar(h5file, "/meta/sys_size",
                        par.L, H5T_NATIVE_DOUBLE);
    }
  }

  // 0 no connection, 1 successfull connection, 2 intersec but no connection
  void write_connectivty_matrix(vec <neuron *> &neurons, hid_t h5file) {
    hid_t h5_matrix = hdf5_create_appendable_nd(
                          h5file, "/connectivity_matrix", H5T_NATIVE_HSIZE,
                          neurons.size(), neurons.size());
    for (size_t i = 0; i < neurons.size(); i++) {
      vec <size_t> row(neurons.size(), 0);
      for (size_t j = 0; j < neurons[i]->outgoing.size(); j++)
        row[ neurons[i]->outgoing[j]->id ] = 1;
      for (size_t j = 0; j < neurons[i]->rejected.size(); j++)
        row[ neurons[i]->rejected[j]->id ] = 2;
      hdf5_append_nd(h5_matrix, row, H5T_NATIVE_HSIZE, i, 0);
    }
    H5Dclose(h5_matrix);
  }

  void write_neuron_details(vec <neuron *> &neurons, hid_t h5file) {
    if (h5file < 0) {
      printf("invalid h5file, no neuron details written");
      return;
    } else {
      vec <real_t> x(neurons.size(), 0);
      vec <real_t> y(neurons.size(), 0);
      vec <real_t> l(neurons.size(), 0);
      vec <real_t> d(neurons.size(), 0);
      for (size_t i = 0; i < neurons.size(); i++) {
        x[i] = neurons[i]->x;
        y[i] = neurons[i]->y;
        l[i] = neurons[i]->l;
        d[i] = neurons[i]->R_d;
      }

      hdf5_write_vector(h5file, "/neurons/pos_x",
                        x, H5T_CUSTOM_REAL);
      hdf5_write_vector(h5file, "/neurons/pos_y",
                        y, H5T_CUSTOM_REAL);
      hdf5_write_vector(h5file, "/neurons/axon_length",
                        l, H5T_CUSTOM_REAL);
      hdf5_write_vector(h5file, "/neurons/radius_dendritic_tree",
                        d, H5T_CUSTOM_REAL);
    }
  }
};

#endif
