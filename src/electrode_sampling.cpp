#ifndef FPS_ELECTRODE_SAMPLING
#define FPS_ELECTRODE_SAMPLING

#include <fstream>
#include <sstream>        // std::stringstream, std::stringbuf
#include <iostream>       // std::cout, std::fixed
#include <iomanip>        // std::setprecision
#include <algorithm>      // std::remove_if
#include <cmath>
#include <vector>
#include <list>
#include <deque>
#include <random>
#include <assert.h>
#include <unistd.h>

#include "rng.h"
#include "helper_io.hpp"
#include "helper_calc.hpp"

// choose the precision for custom float
typedef float real_t;
#define H5T_CUSTOM_REAL (H5T_NATIVE_FLOAT)
// typedef double real_t;
// #define H5T_CUSTOM_REAL (H5T_NATIVE_DOUBLE)

// save some space by aliasing the vector class
template<class T1>
using vec = std::vector <T1>;

// ------------------------------------------------------------------ //
// sampling
// ------------------------------------------------------------------ //

class electrode_sampling {
 public:
  // global parameters that are not part of the topology
  struct parameters {
    // electrode paramters
    size_t NE     = 64;     // number of electrodes
    double d_E    = 100.0;  // [um] dist. between electrodes
    double d_zone = 10.0;   // [um] dist. (radius) of dead zone around elec.
    double gamma  = -1.0;   // dist-scaling of spike contrib. to elec: r^gamma
    size_t cache  = 1e5;    // cache_size before histories are written to disk
  } par;

  vec <electrode *> electrodes;    // vector of all electrodes

  void print_parameters() {
    printf("electrode sampling parameters:\n");
    printf("\tNE       %lu [electrodes] \n", par.NE);
    printf("\td_E      %.2e [um]\n", par.d_E);
    printf("\td_zone   %.2e [um]\n", par.d_zone);
    printf("\tgamma    %.2e\n", par.gamma);
    printf("\tcache    %.0e\n", double(par.cache));
    if (electrodes.size() > 0) {
      printf("\textend   %.2e %.2e | %.2e %.2e [um]\n",
             electrodes.front()->x, electrodes.front()->y,
             electrodes.back()->x,  electrodes.back()->y);
    }
  }

  // signal recording
  size_t at_last = 0;
  vec <size_t> at_history;           // global activity trace
  vec < vec <real_t> > cs_histories;  // coarse time trace
  vec < vec <size_t> > ss_histories;  // subs spike times
  vec < vec <size_t> > at_hist_2d;    // live estimation of m

  // signal writing to hdf5
  hid_t h5dset_cs, h5dset_ss, h5dset_at;
  vec <size_t> cs_offset, ss_offset;

  void init_electrodes(double xy_offset = 0.0) {
    printf("setting up electrode sampling\n");

    // create electrodes
    size_t ne = 0;
    for (size_t i = 0; i < sqrt(par.NE); i++) {
      for (size_t j = 0; j < sqrt(par.NE); j++) {
        electrode *e = new electrode(xy_offset + j*par.d_E,
                                     xy_offset + i*par.d_E, ne);
        electrodes.push_back(e);
        ne += 1;
        // init electrode histories
        ss_histories.push_back(vec <size_t>());
        ss_histories.back().reserve(par.cache);
        cs_histories.push_back(vec <real_t>());
        cs_histories.back().reserve(par.cache);
      }
    }
    at_history.reserve(par.cache);
  }

  // set the neuron contributions based on the topology
  template <typename TOPOLOGY>
  void set_contributions(vec <neuron *> &neurons, TOPOLOGY &tp) {
    printf("setting contributions to electrodes\n");

    for (size_t i = 0; i < neurons.size(); i++) {
      neuron *src = neurons[i];
      // electrode contributions and (alternative) minimum distance fix
      src->electrode_contributions.reserve(electrodes.size());
      for (size_t k = 0; k < electrodes.size(); k++) {
        electrode *e = electrodes[k];
        double rik_squ = tp.get_dist_squ(src, e);
        // this should not happen since we place electrodes first.
        // if neuron too close to electrode, abort
        if (rik_squ < par.d_zone*par.d_zone) {
          printf("neuron too close to electrode.\n");
          exit(-1);
        }

        // coarse sampling contribution
        double irik = pow(rik_squ, par.gamma/2.);
        src->electrode_contributions.push_back(irik);

        // subsampling closest neurons
        if (rik_squ < e->closest_neuron_distance_squ) {
          e->closest_neuron_distance_squ = rik_squ;
          e->closest_neuron = src;
        }
      }
    }
  }

  void init_writing_hdf5(hid_t h5file) {
    if (h5file < 0) {
      h5dset_cs = -1;
      h5dset_ss = -1;
      h5dset_at = -1;
      printf("invalid h5file, failed to setup electrode sampling");
      exit(-1);
    } else {
      // create data sets
      h5dset_cs = hdf5_create_appendable_nd(
                      h5file, "/data/coarse", H5T_CUSTOM_REAL, par.NE,
                      par.cache);
      h5dset_ss = hdf5_create_appendable_nd(
                      h5file, "/data/sub", H5T_NATIVE_HSIZE, par.NE,
                      par.cache);
      h5dset_at = hdf5_create_appendable(
                      h5file, "/data/activity", H5T_NATIVE_HSIZE,
                      par.cache);
      // remember length of each written electrode history (row)
      cs_offset.resize(par.NE, 0);
      ss_offset.resize(par.NE, 0);
    }
  }

  void write_sampling_details(hid_t h5file = -1) {
    if (h5file < 0) {
      printf("invalid h5file, no system details written");
      return;
    } else {
      hdf5_write_scalar(h5file, "/meta/num_elec",
                        par.NE, H5T_NATIVE_HSIZE);
      hdf5_write_scalar(h5file, "/meta/elec_dist",
                        par.d_E, H5T_NATIVE_DOUBLE);
      hdf5_write_scalar(h5file, "/meta/elec_dead_zone",
                        par.d_zone, H5T_NATIVE_DOUBLE);
      hdf5_write_scalar(h5file, "/meta/elec_contribution_exponent",
                        par.gamma, H5T_NATIVE_DOUBLE);

      vec <real_t> x(electrodes.size(), 0);
      vec <real_t> y(electrodes.size(), 0);
      for (size_t i = 0; i < electrodes.size(); i++) {
        x[i] = electrodes[i]->x;
        y[i] = electrodes[i]->y;
      }

      hdf5_write_vector(h5file, "/electrodes/pos_x",
                        x, H5T_CUSTOM_REAL);
      hdf5_write_vector(h5file, "/electrodes/pos_y",
                        y, H5T_CUSTOM_REAL);
    }
  }

  // write cached electrode histories, ideally call this when they match chunks
  void write_histories() {
    for (size_t i = 0; i < electrodes.size(); i++) {
      hdf5_append_nd(h5dset_cs, cs_histories[i], H5T_CUSTOM_REAL, i,
                     cs_offset[i]);
      hdf5_append_nd(h5dset_ss, ss_histories[i], H5T_NATIVE_HSIZE, i,
                     ss_offset[i]);
      cs_offset[i] += cs_histories[i].size();
      ss_offset[i] += ss_histories[i].size();
      cs_histories[i].resize(0);
      ss_histories[i].resize(0);
    }
    hdf5_append(h5dset_at, at_history, H5T_NATIVE_HSIZE);
    at_history.resize(0);
  }

  // also resets electrodes (needed for thermalization)
  inline void measure_step(size_t time = 0, size_t num_active = 0,
                           bool record = true) {
    if (record) {
      update_act_hist(at_hist_2d, at_last, num_active);
      at_history.push_back(num_active);
    }
    for (size_t i = 0; i < electrodes.size(); i++) {
      electrode *e = electrodes[i];
      if (record) {
        // only convert to lower precision float here, when storing,
        // not for individual neuron contributions to electrode signal
        real_t r = real_t(e->cs_signal);
        cs_histories[i].push_back(r);
        if (e->ss_signal) ss_histories[i].push_back(time);
      }
      e->cs_signal = 0.;
      e->ss_signal = false;
    }
    at_last = num_active;
  }

  ~electrode_sampling() {
    if (h5dset_cs > 0) H5Dclose(h5dset_cs);
    if (h5dset_ss > 0) H5Dclose(h5dset_ss);
    if (h5dset_at > 0) H5Dclose(h5dset_at);
    for (size_t i = 0; i < electrodes.size(); i++) delete electrodes[i];
  }
};

// ------------------------------------------------------------------ //
// dynamics and updates of electrodes
// ------------------------------------------------------------------ //

class dynamic_branching {
 public:
  // global parameters
  struct parameters {
    double m      = 0.98;   // microscopic branching param., also in topology
    double h      = 4e-5;   // drive (probability for spontaneous activation)
  } par;

  // we need the neuron (and electrode) vector in order to update things
  vec <neuron *>    &neurons;
  vec <electrode *> &electrodes;

  // keep track of active neurons
  std::deque <neuron *> active_neurons;
  size_t num_active_new = 0;
  size_t num_active_old = 0;

  dynamic_branching(vec <neuron *> &n_,
                    vec <electrode *> &e_, double m_, double h_) :
    neurons(n_), electrodes(e_) {
    par.m = m_;
    par.h = h_;
    printf("setting up dynamics\n\tcoalesence compensation: ");
    #ifdef NOCC
    printf("No\n");
    #else
    printf("Yes\n");
    #endif
    printf("\tm: %.3e\n", par.m);
    printf("\th: %.3e\n", par.h);
  }

  // multiply all outgoing probabilities with par.m
  void scale_activation_probabilities(vec <neuron *> &neurons) {
    printf("scaling activation probabilities to m = %.2f\n", par.m);
    for (size_t i = 0; i < neurons.size(); i++) {
      neuron *src = neurons[i];
      for (size_t j = 0; j < src->outgoing_probability.size(); j++) {
        src->outgoing_probability[j] *= par.m;
      }
    }

    // consistency checks and some insight into how local connectivity is
    // integrated probability, should equal one
    double k0 = 0., k1=0., k2=0., k3=0., p = 0.;
    for (size_t i = 0; i < neurons.size(); i++) {
      neuron *src = neurons[i];
      double p_temp = 0.;
      for (size_t j = 0; j < src->outgoing.size(); j++) {
        if (src->outgoing_probability.size() == 1) {
          // all connections share same probability (orlandi topology)
          p_temp += src->outgoing_probability[0];
          // printf("%lu %lu %.4f\n", i, j, src->outgoing_probability[0]);
        } else {
          // different prob for every connection (local gauss topology)
          p_temp += src->outgoing_probability[j];
          // printf("%lu %lu %.4f\n", i, j, src->outgoing_probability[0]);
        }
        if (p_temp < 1./3./par.m) k0 += 1;
        if (p_temp < 1./2./par.m) k1 += 1;
        if (p_temp < 2./3./par.m) k2 += 1;
        k3 += 1;
      }
      p += p_temp;
    }
    p  /= double(neurons.size());
    k0 /= double(neurons.size());
    k1 /= double(neurons.size());
    k2 /= double(neurons.size());
    k3 /= double(neurons.size());
    printf("\tintegrated probability per neuron: %.2f\n", p);
    printf("\tnumber of connections responsible for probability:\n");
    printf("\t\tp 1/3 ~ %.2f | p 1/2 ~ %.2f | p 2/3 ~ %.2f | p 1 ~ %.2f \n",
           k0, k1, k2, k3);
    // would be good to do the same for the distance to the kth neuron
  }

  inline void activate_neuron(neuron *n) {
    // manually check not to activate already active neurons!
    n->active = true;
    num_active_new += 1;
    active_neurons.push_back(n);
    // update electrodes
    for (size_t i = 0; i < electrodes.size(); i++) {
      electrode *e = electrodes[i];
      e->cs_signal += n->electrode_contributions[i];
      if (e->closest_neuron == n) e->ss_signal = true;
    }
  }

  inline void update_step() {
    size_t temp = 0;
    update_step(temp);
  }

  inline void update_step(size_t &num_active) {

    num_active_old     = num_active_new;
    num_active_new     = 0;

    // spontanious activation and resetting
    for (size_t i = 0; i < neurons.size(); i++) {
      if (udis(rng) < par.h) {
        activate_neuron(neurons[i]);
      } else {
        neurons[i]->active = false;
      }
    }

    // spread activity. if active, no activation possible
    for (size_t i = 0; i < num_active_old; i++) {
      neuron *src = active_neurons[i];
      size_t num_recurrent = 0;
      for (size_t j = 0; j < src->outgoing.size(); j++) {
        double p_threshold;
        if (src->outgoing_probability.size() == 1)
          // same probability for all outgoing connections
          p_threshold = src->outgoing_probability[0];
        else
          // different probability for every outgoing connection
          p_threshold = src->outgoing_probability[j];
        #ifdef NOCC
        // simple recurrent activations. if target already active, skip.
        if (!src->outgoing[j]->active
            && udis(rng) < p_threshold) {
          activate_neuron(src->outgoing[j]);
        }
        #else
        // coalesence compensation
        // if random number says to activate, try current target.
        // if current target already active, try the next closest (farther)
        // connected neuron until successful or no more candidates
        if (udis(rng) < p_threshold) {
          if (!src->outgoing[j]->active) {
            activate_neuron(src->outgoing[j]);
            num_recurrent += 1;
          } else if (num_recurrent < src->outgoing.size()) {
            bool reassigned = false;
            while (!reassigned && j+1 < src->outgoing.size()) {
              j += 1;
              if (!src->outgoing[j]->active) {
                activate_neuron(src->outgoing[j]);
                num_recurrent += 1;
                reassigned = true;
              }
            }
          } else {
            printf("\n\n num recurrenct reached limit! src id: %lu\n\n",
                   src->id);
          }
        }
        #endif

      }
    }
    assert(active_neurons.size() == num_active_new+num_active_old);

    // deactivate by cleaning list of active neurons
    active_neurons.erase(active_neurons.begin(),
                         active_neurons.begin() + num_active_old);

    assert(active_neurons.size() == num_active_new);
    num_active = num_active_new;
  }

  void write_dynamic_details(hid_t h5file) {
    if (h5file < 0) {
      printf("invalid h5file, no dynamic details written");
      return;
    } else {
      #ifdef NOCC
      size_t cc = 0;
      #else
      size_t cc = 1;
      #endif
      hdf5_write_scalar(h5file, "/meta/coalesence_compensation_bool",
                        cc, H5T_NATIVE_HSIZE);
      hdf5_write_scalar(h5file, "/meta/m_micro",
                        par.m, H5T_NATIVE_DOUBLE);
      hdf5_write_scalar(h5file, "/meta/h_prob",
                        par.h, H5T_NATIVE_DOUBLE);
    }
  }
};
#endif
