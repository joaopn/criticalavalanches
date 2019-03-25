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

#include "helper_io.hpp"
#include "helper_calc.hpp"

// rng and distributions
std::mt19937 rng(316);
std::uniform_real_distribution<> udis(0, 1);
std::normal_distribution<> gdis(0, 1);
std::binomial_distribution<int> bdis(10, 0.5);

class electrode;
class neuron {
 public:
  bool   active;
  double x;
  double y;
  size_t id;
  std::vector< neuron * >    outgoing;
  std::vector< electrode * > electrodes;
  std::vector< double >     outgoing_probability;
  std::vector< double >     electrode_contributions;

  neuron(double x_, double y_, double id_ = 0) {
    active = false;
    x = x_;
    y = y_;
    id = id_;
  }
};

class electrode {
 public:
  double cs_signal;
  bool   ss_signal;
  double x;
  double y;
  size_t id;
  neuron *closest_neuron;
  double closest_neuron_distance_squ;

  electrode(double x_, double y_, double id_ = 0) {
    x = x_;
    y = y_;
    id = id_;
    closest_neuron_distance_squ = std::numeric_limits<double>::max();
    cs_signal = 0.;
    ss_signal = false;
  }
};

class asystem {
 public:

  // ------------------------------------------------------------------ //
  // system variables
  // ------------------------------------------------------------------ //
  size_t num_neur;
  size_t num_outgoing;
  size_t num_elec;
  double sys_size;
  double elec_frac;
  double delta_l;
  double delta_l_nn;
  double neuron_density;
  double outgoing_distance;
  double m_micro;
  double h_prob;
  double w_var;
  std::vector< neuron * >    neurons;
  std::vector< electrode * > electrodes;

  // keep track of active neurons
  std::deque < neuron * > active_neurons;
  size_t num_active_new;
  size_t num_active_old;

  // signal recording
  size_t cache_size;
  std::vector< size_t > at_history;                   // global activity trace
  std::vector< std::vector< double > > cs_histories;  // coarse time trace
  std::vector< std::vector< size_t > > ss_histories;  // subs spike times
  std::vector< std::vector< size_t > > at_hist_2d;    // live estimation of m

  // ------------------------------------------------------------------ //
  // construct system
  // ------------------------------------------------------------------ //
  asystem(double sys_size_ = 1., double elec_frac_ = .25,
          size_t num_neur_ = 144000, size_t num_outgoing_ = 1000,
          size_t num_elec_ = 100, double m_micro_ = 1.0,
          double w_var_ = 1., double h_prob_ = 0.1) {
    num_neur       = num_neur_;
    num_outgoing   = num_outgoing_;
    num_elec       = num_elec_;
    sys_size       = sys_size_;
    elec_frac      = elec_frac_;
    m_micro        = m_micro_;
    w_var          = w_var_;
    h_prob         = h_prob_;
    neuron_density = double(num_neur)/sys_size/sys_size;
    num_active_old = 0;
    num_active_new = 0;
    cache_size     = size_t(1e4);

    // analytic solution for average inter neuron distance delta_l
    delta_l = pow(sys_size, 3.)/6. * ( sqrt(2.) + log(1. + sqrt(2.)) );

    // analytic solution for average nearest-neighbour distance
    delta_l_nn = -exp(-neuron_density * M_PI)
                 + erf(sqrt(neuron_density * M_PI))/2./sqrt(neuron_density);

    #ifndef COALCOMP
    printf("creating system\n");
    #else
    printf("creating system with coalesence compensation\n");
    #endif
    printf("\tnumber of neurons: %lu\n", num_neur);
    printf("\th: %.3e\n", h_prob);
    printf("\tm: %.3e\n", m_micro);
    printf("\tw: %.3e\n", w_var);

    // create neurons
    for (size_t i = 0; i < num_neur; i++) {
      neuron *n = new neuron(udis(rng)*sys_size, udis(rng)*sys_size, i);
      neurons.push_back(n);
    }

    // choose r so that approx num_outgoing neighbours are accessible
    outgoing_distance = sqrt(sys_size*sys_size*num_outgoing/M_PI
                             /double(num_neur));

    // create connections
    double mc = connect_neurons_using_interaction_radius(outgoing_distance);

    printf("\toutgoing connections per neuron: ~%lu\n", num_outgoing);
    #ifndef NDEBUG
    printf("\t\t(measured: %.2f)\n", mc);
    #endif
    printf("\tconnection distance: %.2e\n", outgoing_distance);
    printf("\taverage distance between neurons: %.2e\n", delta_l);
    #ifndef NDEBUG
    printf("\t\t(measured: %.2e)\n", measure_avg_distance());
    #endif
    printf("\taverage distance between nearest neighbours: %.2e\n",
           delta_l_nn);
    #ifndef NDEBUG
    printf("\t\t(measured: %.2e)\n", measure_avg_nearest_neighbour_distance());
    #endif

    // create electrodes to spread over a frac of the system (each direction)
    size_t ne = 0;
    double de = sys_size*elec_frac/sqrt(double(num_elec));
    for (size_t i = 0; i < sqrt(num_elec); i++) {
      for (size_t j = 0; j < sqrt(num_elec); j++) {
        electrode *e = new electrode(j*de, i*de, ne);
        electrodes.push_back(e);
        ne += 1;
        // init electrode histories
        ss_histories.push_back(std::vector< size_t >());
        ss_histories.back().reserve(cache_size);
        cs_histories.push_back(std::vector< double >());
        cs_histories.back().reserve(cache_size);
      }
    }
    at_history.reserve(cache_size);

    printf("electrodes placed\n");
    printf("\tnumber of electrodes: %lu^2 = %lu\n",
           size_t(sqrt(num_elec)), ne);
    printf("\telectrode distance: %.2e\n", de);

    // precalculate how much each neuron's spike contributes to an electrode
    // and set the probabilities for recurrent activations
    set_contributions_and_probabilities(de);
  }

  // ------------------------------------------------------------------ //
  // topology
  // ------------------------------------------------------------------ //
  // find distance between neurons (or electrodes)
  template <typename T1, typename T2>
  inline double get_dist_squ(T1 *p, T2 *t) {
    double xr = std::max(t->x,p->x);
    double xl = std::min(t->x,p->x);
    double yr = std::max(t->y,p->y);
    double yl = std::min(t->y,p->y);

    double x2 = std::min((xr-xl)*(xr-xl), (sys_size-xr+xl)*(sys_size-xr+xl));
    double y2 = std::min((yr-yl)*(yr-yl), (sys_size-yr+yl)*(sys_size-yr+yl));

    return (x2+y2);
  }

  // create connections with hard limit given by radius
  double connect_neurons_using_interaction_radius(double radius) {
    printf("connecting neurons with radius %.2e\n", radius);
    double dist_squ_limit = radius*radius;
    size_t avg_connection_count = 0;
    for (size_t i = 0; i < neurons.size(); i++) {
      neuron *src = neurons[i];
      src->outgoing.reserve(size_t(num_outgoing*1.1));
      for (size_t j = 0; j < neurons.size(); j++) {
        neuron *tar = neurons[j];
        if (src == tar) continue;
        if (get_dist_squ(src, tar) < dist_squ_limit) {
          src->outgoing.push_back(tar);
          avg_connection_count += 1;
        }
      }
      if(i==0 || is_percent(i, size_t(neurons.size()), 10.)) {
        printf("\t%s, %lu/%lu\n", time_now().c_str(), i, neurons.size());
      }
    }
    printf("done\n");
    return avg_connection_count/double(neurons.size());
  }

  // effective interaction radius is set in here via w_var
  void set_contributions_and_probabilities(double elec_distance) {
    printf("calculating interaction radius and contributions to electrodes\n");
    printf("\tw_eff ~%.1e r_max\n", w_var*delta_l_nn/outgoing_distance);
    // calculate contributions to electrodes and probabilities to activate
    double w_squ       = 2.*pow(w_var*delta_l_nn, 2.);
    double dik_min_squ =    pow(elec_distance/10., 2.);
    // if we do not ensure this, neurons will always be too close to some elec
    assert(dik_min_squ < pow(.5*elec_distance, 2.));

    for (size_t i = 0; i < neurons.size(); i++) {
      if(i==0 || is_percent(i, size_t(neurons.size()), 10.)) {
        printf("\t%s, %lu/%lu\n", time_now().c_str(), i, neurons.size());
      }
      neuron *src = neurons[i];

      // electrode contributions and minimum distance fix
      src->electrode_contributions.reserve(electrodes.size());
      for (size_t k = 0; k < electrodes.size(); k++) {
        electrode *e = electrodes[k];
        double dik_squ = get_dist_squ(src, e);
        // if neuron too close to electrode, move neuron a bit
        if (dik_squ < dik_min_squ) {
          double dx = src->x - e->x;
          double dy = src->y - e->y;
          double a = dy/dx;
          double r = 1.01*pow(dik_min_squ/dik_squ, .5);
          double xn = e->x + dx*r;
          double yn = e->y + a*dx*r;
          src->x = xn-std::floor(xn/sys_size)*sys_size;
          src->y = yn-std::floor(yn/sys_size)*sys_size;
          dik_squ = get_dist_squ(src, e);
          assert(dik_squ > dik_min_squ);
        }

        double idik = pow(dik_squ, -1./2.);
        src->electrode_contributions.push_back(idik);
        // subsampling closest neurons
        if (dik_squ < e->closest_neuron_distance_squ) {
          e->closest_neuron_distance_squ = dik_squ;
          e->closest_neuron = src;
        }
      }

      // activation probabilities based on distance
      size_t n_cout = src->outgoing.size();
      src->outgoing_probability.reserve(n_cout);
      double norm = 0.;
      for (size_t j = 0; j < n_cout; j++) {
        double dij_squ = get_dist_squ(src, src->outgoing[j]);
        double pij = exp(-dij_squ/w_squ);
        src->outgoing_probability.push_back(pij);
        norm += pij;
      }

      // normalize probabilities to m_micro
      for (size_t j = 0; j < n_cout; j++) {
        src->outgoing_probability[j] *= m_micro/norm;
      }
    }
    printf("done\n");
  }

  // delta_l_nn
  double measure_avg_nearest_neighbour_distance() {
    double avg_dist = 0;
    for (size_t i = 0; i < neurons.size(); i++) {
      double shortest_distance = std::numeric_limits<double>::max();
      for (size_t j = 0; j < neurons[i]->outgoing.size(); j++) {
        // if (neurons[i] == j) continue;
        assert(neurons[i] != neurons[i]->outgoing[j]);
        double d = get_dist_squ(neurons[i], neurons[i]->outgoing[j]);
        if (d < shortest_distance) shortest_distance = d;
      }
      avg_dist += std::sqrt(shortest_distance);
    }
    avg_dist = avg_dist / double(neurons.size());
    return avg_dist;
  }

  double measure_avg_connection_length() {
    double avg_dist = 0;
    for (size_t i = 0; i < neurons.size(); i++) {
      for (size_t j = 0; j < neurons[i]->outgoing.size(); j++) {
        double d = get_dist_squ(neurons[i], neurons[i]->outgoing[j]);
        avg_dist += std::sqrt(d);
      }
    }
    avg_dist = avg_dist / double(num_neur) / double(num_outgoing);
    return avg_dist;
  }

  // delta_l
  double measure_avg_distance() {
    double dl = 0;
    for (size_t i = 0; i < neurons.size(); i++) {
      for (size_t j = 0; j < neurons.size(); j++) {
        double d = get_dist_squ(neurons[i], neurons[j]);
        dl += std::sqrt(d);
      }
    }
    dl = dl / double(neurons.size()) / double(neurons.size());
    return dl;
  }

  // ------------------------------------------------------------------ //
  // simulation
  // ------------------------------------------------------------------ //
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

    num_active_old = num_active_new;
    num_active_new = 0;

    // spontanious activation and resetting
    for (size_t i = 0; i < neurons.size(); i++) {
      if (udis(rng) < h_prob) {
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
        #ifndef COALCOMP
        // default recurrent activations. if target already active, skip.
        if (!src->outgoing[j]->active
            && udis(rng) < src->outgoing_probability[j]) {
          activate_neuron(src->outgoing[j]);
        }
        #else
        // coalesence compensation, naively implemented.
        // if random number says to activate, try current target.
        // if current target already active, randomly try other connected
        // neurons until successful
        // this breaks locality a bit
        if (udis(rng) < src->outgoing_probability[j]) {
          if (!src->outgoing[j]->active) {
            activate_neuron(src->outgoing[j]);
            num_recurrent += 1;
          } else if (num_recurrent < src->outgoing.size()) {
            bool reassigned = false;
            while (!reassigned) {
              size_t new_j = floor(udis(rng)*src->outgoing.size());
              if (!src->outgoing[new_j]->active) {
                activate_neuron(src->outgoing[new_j]);
                num_recurrent += 1;
                reassigned = true;
              }
            }
          } else {
            printf("avoid this!\n");
            exit(-1);
          }
        }
        #endif

      }
    }
    assert(active_neurons.size() == num_active_new+num_active_old);

    // deactivate by cleaning list of active neurons
    // for (size_t i = 0; i < num_active_old; i++) {
    //   active_neurons.pop_front();
    // }
    active_neurons.erase(active_neurons.begin(),
                         active_neurons.begin() + num_active_old);

    assert(active_neurons.size() == num_active_new);
  }

  inline void measure_step(size_t time = 0, bool record = true) {
    if (record) {
      at_history.push_back(num_active_new);
      update_act_hist(at_hist_2d, num_active_old, num_active_new);
    }
    for (size_t i = 0; i < electrodes.size(); i++) {
      electrode *e = electrodes[i];
      if (record) {
        cs_histories[i].push_back(e->cs_signal);
        if (e->ss_signal) ss_histories[i].push_back(time);
      }
      e->cs_signal = 0.;
      e->ss_signal = false;
    }
  }
};


class exporter {
 public:
  asystem *sys;
  std::string h5filepath, h5dset_csname, h5dset_ssname;
  hid_t h5file, h5dset_cs, h5dset_ss, h5dset_at;
  std::vector< size_t > cs_offset, ss_offset;

  exporter(std::string filepath, asystem *sys_, size_t seed) {
    sys = sys_;

    if (filepath.rfind('/')!=std::string::npos) {
      std::string dirname = filepath;
      dirname.erase(dirname.rfind('/'));
      system(("mkdir -p " + dirname).c_str());
    }
    printf("exporting files to %s\n", filepath.c_str());

    // open file and create groups
    h5file = H5Fcreate(filepath.c_str(), H5F_ACC_TRUNC,
                       H5P_DEFAULT, H5P_DEFAULT);

    if(h5file < 0) exit(-1);

    H5Gcreate(h5file, "/meta", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gcreate(h5file, "/data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // write system details
    hdf5_write_scalar(h5file, "/meta/seed",
                      seed, H5T_NATIVE_HSIZE);
    hdf5_write_scalar(h5file, "/meta/num_neur",
                      sys->num_neur, H5T_NATIVE_HSIZE);
    hdf5_write_scalar(h5file, "/meta/num_outgoing",
                      sys->num_outgoing, H5T_NATIVE_HSIZE);
    hdf5_write_scalar(h5file, "/meta/num_elec",
                      sys->num_elec, H5T_NATIVE_HSIZE);
    hdf5_write_scalar(h5file, "/meta/sys_size",
                      sys->sys_size, H5T_NATIVE_DOUBLE);
    hdf5_write_scalar(h5file, "/meta/elec_frac",
                      sys->elec_frac, H5T_NATIVE_DOUBLE);
    hdf5_write_scalar(h5file, "/meta/delta_l",
                      sys->delta_l, H5T_NATIVE_DOUBLE);
    hdf5_write_scalar(h5file, "/meta/delta_l_nn",
                      sys->delta_l_nn, H5T_NATIVE_DOUBLE);
    hdf5_write_scalar(h5file, "/meta/neuron_density",
                      sys->neuron_density, H5T_NATIVE_DOUBLE);
    hdf5_write_scalar(h5file, "/meta/outgoing_distance",
                      sys->outgoing_distance, H5T_NATIVE_DOUBLE);
    hdf5_write_scalar(h5file, "/meta/w_var",
                      sys->w_var, H5T_NATIVE_DOUBLE);
    hdf5_write_scalar(h5file, "/meta/m_micro",
                      sys->m_micro, H5T_NATIVE_DOUBLE);
    hdf5_write_scalar(h5file, "/meta/h_prob",
                      sys->h_prob, H5T_NATIVE_DOUBLE);

    // create data sets
    h5dset_cs = hdf5_create_appendable_nd(
                    h5file, "/data/coarse", H5T_NATIVE_DOUBLE, sys->num_elec,
                    sys->cache_size);
    h5dset_ss = hdf5_create_appendable_nd(
                    h5file, "/data/sub", H5T_NATIVE_HSIZE, sys->num_elec,
                    sys->cache_size);
    h5dset_at = hdf5_create_appendable(
                    h5file, "/data/activity", H5T_NATIVE_HSIZE,
                    sys->cache_size);

    // remember length of each written electrode history (row)
    cs_offset.resize(sys->num_elec, 0);
    ss_offset.resize(sys->num_elec, 0);
  }

  void finalize() {
    H5Fclose(h5file);
    H5Dclose(h5dset_cs);
    H5Dclose(h5dset_ss);
    H5Dclose(h5dset_at);
  }

  // write cached electrode histories, ideally call this when they match chunks
  void write_histories() {
    for (size_t i = 0; i < sys->electrodes.size(); i++) {
      hdf5_append_nd(h5dset_cs, sys->cs_histories[i], H5T_NATIVE_DOUBLE, i,
                     cs_offset[i]);
      hdf5_append_nd(h5dset_ss, sys->ss_histories[i], H5T_NATIVE_HSIZE, i,
                     ss_offset[i]);
      cs_offset[i] += sys->cs_histories[i].size();
      ss_offset[i] += sys->ss_histories[i].size();
      sys->cs_histories[i].resize(0);
      sys->ss_histories[i].resize(0);
    }
    hdf5_append(h5dset_at, sys->at_history, H5T_NATIVE_HSIZE);
    sys->at_history.resize(0);
  }

  // write system configuration. vectors with ptrs to electrodes or neurons
  template <typename T1>
  void write_config(std::string filename, std::vector<T1> &config_to_save) {
    std::ofstream file (filename, std::ofstream::out);
    file << std::setprecision(0) << std::fixed;
    file << "#L=" << sys->sys_size <<std::endl;
    file << "#N_exported=" << config_to_save.size() <<std::endl;
    file << "#x y" << std::endl;
    file << std::setprecision(7) << std::scientific;
    for (size_t i = 0; i<config_to_save.size(); i++) {
      file << config_to_save[i]->x << "\t" << config_to_save[i]->y;
      file << std::endl;
    }
    file.close();
    file.clear();
  }
};

// ------------------------------------------------------------------ //
// main
// ------------------------------------------------------------------ //
int main(int argc, char *argv[]) {

  double sys_size     = 1.;         // sytem size
  double time_steps   = 1e3;        // number of time steps
  size_t num_neur     = 144000;     // number of neurons
  size_t num_outgoing = 1000;       // average outgoing connections per neuron
  size_t num_elec     = 100;        // total number of electrodes
  double elec_frac    = .25;        // electrodes cover part of system per dim
  size_t seed         = 314;        // seed for the random number generator
  double m_micro      = 1.0;        // branching parameter applied locally
  double w            = 4.;         // eff. conn-length [unit=nearestneur-dist]
  double h            = .01;        // probability for spontaneous activation
  std::string path    = "";         // output path for results

  for (size_t i=0; i < argc; i++) {
    if (i+1 != argc) {
      if(std::string(argv[i])=="-T") time_steps     = atof(argv[i+1]);
      if(std::string(argv[i])=="-N") num_neur       = atof(argv[i+1]);
      if(std::string(argv[i])=="-k") num_outgoing   = atof(argv[i+1]);
      if(std::string(argv[i])=="-e") num_elec       = atof(argv[i+1]);
      if(std::string(argv[i])=="-f") elec_frac      = atof(argv[i+1]);
      if(std::string(argv[i])=="-s") seed           = atof(argv[i+1]);
      if(std::string(argv[i])=="-m") m_micro        = atof(argv[i+1]);
      if(std::string(argv[i])=="-w") w              = atof(argv[i+1]);
      if(std::string(argv[i])=="-h") h              = atof(argv[i+1]);
      if(std::string(argv[i])=="-o") path           =      argv[i+1] ;
    }
  }

  if (path == "") {
    printf("specify output path with '-o'\n");
    return -1;
  }

  setbuf(stdout, NULL); // print direct without calling flush, also on cluster
  rng.seed(1000+seed);
  rng.discard(70000);

  asystem *sys = new asystem(sys_size, elec_frac, num_neur, num_outgoing,
                             num_elec, m_micro, w, h);
  exporter *out = new exporter(path, sys, seed);

  // thermalization with ~1/h time steps, dont run with h=0
  // printf("thermalizing for %.0e time steps\n", 1./h);
  // for (size_t i = 0; i < size_t(1./h); i++) {
  //   if(is_percent(i, size_t(1./h), 10.) || have_passed_hours(6.)) {
  //     printf("\t%s, progress ~%2.0f%%, activity ~%.3f\n",
  //            time_now().c_str(), is_percent(i, size_t(1./h), 10.),
  //            sys->num_active_old/double(sys->num_neur));
  //   }
  //   sys->update_step();
  //   sys->measure_step(i, false);
  // }


  printf("simulating for %.0e time steps\n", time_steps);
  for (size_t i = 0; i < size_t(time_steps); i++) {
    if(is_percent(i, size_t(time_steps), 10.) || have_passed_hours(6.)) {
      double act, var, mlr;
      var = variance(sys->at_history, act);
      mlr = m_from_lin_regr(sys->at_hist_2d);
      printf("\t%s, ~%2.0f%%, last %.0e steps: act ~%.3f, cv ~%.2f, mlr ~%.5f, tau(2ms) ~%.1f\n",
             time_now().c_str(), is_percent(i, size_t(time_steps), 10.),
             double(sys->at_history.size()),
             // tau assuming timesteps of 2ms
             act/sys->num_neur, sqrt(var)/act, mlr, -2./log(mlr));
      // out->write_histories();
    }
    sys->update_step();
    sys->measure_step(i, true);
    if(i%sys->cache_size==0) out->write_histories();
  }
  out->write_histories();
  printf("done\n\n");
  out->finalize();

  return 0;


  // ------------------------------------------------------------------ //
  // topology test
  // ------------------------------------------------------------------ //
  out->write_config("/Users/paul/Desktop/neurons_2.dat", sys->neurons);
  out->write_config("/Users/paul/Desktop/electrodes_2.dat", sys->electrodes);

  std::vector< neuron * > test(100, nullptr);
  for (size_t i = 0; i < sys->electrodes.size(); i++) {
    electrode *e = sys->electrodes[i];
    printf("%lu %e %e\n", e->id, e->x, e->y);
    neuron *n = e->closest_neuron;
    printf("\t%lu %e %e\n", n->id, n->x, n->y);
    test[i] = n;
  }
  out->write_config("/Users/paul/Desktop/closest_2.dat", test);

  neuron *t = sys->neurons[size_t(udis(rng)*num_neur)];
  std::vector< neuron * > foo = {t};

  double d = 0.;
  for (size_t i = 0; i < t->outgoing.size(); i++) {
    d += t->outgoing_probability[i];
    printf("prob: %e %e\n", t->outgoing_probability[i], d);
  }

  for (size_t i = 0; i < t->electrode_contributions.size(); i++) {
    printf("contrib: %e\n", t->electrode_contributions[i]);
  }

  out->write_config("/Users/paul/Desktop/conn_2.dat", foo);
  out->write_config("/Users/paul/Desktop/conns_2.dat", t->outgoing);

  delete sys;

}
