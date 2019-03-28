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

#include "helper_io.hpp"
#include "helper_calc.hpp"

// choose the precision for time series of electrode signal
typedef float real_t;
#define H5T_CUSTOM_REAL (H5T_NATIVE_FLOAT)
// typedef double real_t;
// #define H5T_CUSTOM_REAL (H5T_NATIVE_DOUBLE)



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
  size_t num_neur;                // number of neurons
  double neur_dist;               // avg. closest distance between neurons
  double sys_size;                // linear system size, derived from above
  size_t num_elec;                // total number of electrodes
  double elec_dist;               // distance between electrodes [neur_dist]
  size_t num_outgoing;            // number outgoing connections per neuron
  double d_max;                   // maximum distance for neur. connections
  double r_min;                   // minimum distance betw. neur. and elec.
  double neuron_density;          // neuron density
  double m_micro;                 // microscopic branching parameter
  double h_prob;                  // probability of spontanious activation
  double sigma;                   // effective connection radius (via gaussian)
  double delta_l;                 // dist between all neurons, system wide

  std::vector< neuron * >    neurons;       // vector of all neurons
  std::vector< electrode * > electrodes;    // vector of all electrodes

  // keep track of active neurons
  std::deque < neuron * > active_neurons;
  size_t num_active_new;
  size_t num_active_old;

  // signal recording
  size_t cache_size;
  std::vector< size_t > at_history;                   // global activity trace
  std::vector< std::vector< real_t > > cs_histories;  // coarse time trace
  std::vector< std::vector< size_t > > ss_histories;  // subs spike times
  std::vector< std::vector< size_t > > at_hist_2d;    // live estimation of m

  // ------------------------------------------------------------------ //
  // construct system
  // ------------------------------------------------------------------ //
  asystem(
      size_t num_neur_, double neur_dist_,
      size_t num_elec_, double elec_dist_,
      size_t num_outgoing_,
      double m_micro_, double sigma_, double h_prob_,
      size_t cache_size_) {
    num_neur       = num_neur_;
    neur_dist      = neur_dist_;
    num_elec       = num_elec_;
    elec_dist      = elec_dist_;
    num_outgoing   = num_outgoing_;

    // find system size mathcing neur_distance numerically
    double delta_l_nn;
    sys_size = 2*sqrt(num_neur)*neur_dist;

    neuron_density = double(num_neur)/sys_size/sys_size;
    // analytic solution for average nearest-neighbour distance
    // sys_size correct if this matches neur_dist
    delta_l_nn = -exp(-neuron_density * M_PI)
                 + erf(sqrt(neuron_density * M_PI))/2./sqrt(neuron_density);

    m_micro        = m_micro_;
    sigma          = sigma_;
    h_prob         = h_prob_;
    cache_size     = cache_size_;

    num_active_old = 0;
    num_active_new = 0;


    // analytic solution for average inter neuron distance delta_l
    delta_l = pow(sys_size, 3.)/6. * ( sqrt(2.) + log(1. + sqrt(2.)) );

    printf("creating system\n\tcoalesence compensation: ");
    #ifdef NOCC
    printf("No\n");
    #else
    printf("Yes\n");
    #endif
    printf("\tnumber of neurons: %lu\n", num_neur);
    printf("\tnearest neuron distance: %.3e\n", neur_dist);
    printf("\tsystem size: %.3e\n", sys_size);
    printf("\th: %.3e\n", h_prob);
    printf("\tm: %.3e\n", m_micro);
    printf("\tsigma: %.3e\n", sigma);


    // create electrodes
    size_t ne = 0;
    double de = elec_dist*neur_dist;
    for (size_t i = 0; i < sqrt(num_elec); i++) {
      for (size_t j = 0; j < sqrt(num_elec); j++) {
        electrode *e = new electrode(j*de, i*de, ne);
        electrodes.push_back(e);
        ne += 1;
        // init electrode histories
        ss_histories.push_back(std::vector< size_t >());
        ss_histories.back().reserve(cache_size);
        cs_histories.push_back(std::vector< real_t >());
        cs_histories.back().reserve(cache_size);
      }
    }
    at_history.reserve(cache_size);

    // dead zones around electrodes where neurons are forbidden
    r_min = neur_dist/5.;
    size_t rejections = 0;

    printf("electrodes placed\n");
    printf("\tnumber of electrodes: %lu^2 = %lu\n",
           size_t(sqrt(num_elec)), ne);
    printf("\telectrode distance: %.2e\n", de);
    printf("\telectrode dead zone (r_min): %.2e\n", r_min);
    printf("\telectrode array size: %.2e\n", elec_dist*sqrt(num_elec));

    // create neurons
    while (neurons.size() < num_neur) {
      bool placed = false;
      neuron *n;
      while (!placed) {
        n = new neuron(udis(rng)*sys_size, udis(rng)*sys_size, neurons.size());
        placed = true;
        for (size_t k = 0; k < electrodes.size(); k++) {
          if (r_min*r_min >= get_dist_squ(electrodes[k], n)) {
            delete n;
            placed = false;
            rejections += 1;
            if (rejections > num_neur/10) {
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

    printf("neurons placed\n\tnumber of neurons: %lu\n", neurons.size());
    printf("\trejections in dead zones: %lu\n", rejections);

    // create connections
    // choose d so that approx num_outgoing neighbours are accessible
    d_max = sqrt(sys_size*sys_size*num_outgoing/M_PI/double(num_neur));
    double mc = connect_neurons_using_interaction_radius(d_max);

    printf("\toutgoing connections per neuron: ~%lu\n", num_outgoing);
    #ifndef NDEBUG
    printf("\t\t(measured: %.2f)\n", mc);
    #endif
    printf("\tconnection distance: %.2e\n", d_max);
    printf("\taverage distance between neurons: %.2e\n", delta_l);
    #ifndef NDEBUG
    printf("\t\t(measured: %.2e)\n", measure_avg_distance());
    #endif
    printf("\taverage distance between nearest neighbours: %.2e\n",
           neur_dist);
    // #ifndef NDEBUG
    // sys_size correct if this matches neur_dist
    printf("\t\t(measured: %e)\n", measure_avg_nearest_neighbour_distance());
    // #endif

    // precalculate how much each neuron's spike contributes to an electrode
    // and set the probabilities for recurrent activations
    set_contributions_and_probabilities();
  }

  ~asystem() {
    for (size_t i = 0; i < neurons.size(); i++) delete neurons[i];
    for (size_t i = 0; i < electrodes.size(); i++) delete electrodes[i];
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

      std::vector< std::pair <double, neuron *> > sortable;
      sortable.reserve(size_t(num_outgoing*1.1));

      for (size_t j = 0; j < neurons.size(); j++) {
        neuron *tar = neurons[j];
        if (src == tar) continue;
        double dist_squ = get_dist_squ(src, tar);
        if (dist_squ < dist_squ_limit) {
          sortable.push_back( std::make_pair(dist_squ, tar) );
          avg_connection_count += 1;
        }
      }

      // sort outgoing connections by distance
      sort(sortable.begin(), sortable.end());
      src->outgoing.reserve(sortable.size());
      for (size_t j = 0; j < sortable.size(); j++) {
        src->outgoing.push_back(sortable[j].second);
      }


      if(i==0 || is_percent(i, size_t(neurons.size()), 10.)) {
        printf("\t%s, %lu/%lu\n", time_now().c_str(), i, neurons.size());
      }
    }
    printf("done\n");
    return avg_connection_count/double(neurons.size());
  }

  // effective interaction radius is set in here via sigma
  void set_contributions_and_probabilities() {
    printf("calculating interaction radius and contributions to electrodes\n");
    printf("\tsigma_eff: %.3e\n", sigma*neur_dist);
    // calculate contributions to electrodes and probabilities to activate
    double sigma_squ   = 2.*pow(sigma*neur_dist, 2.);

    for (size_t i = 0; i < neurons.size(); i++) {
      if(i==0 || is_percent(i, size_t(neurons.size()), 10.)) {
        printf("\t%s, %lu/%lu\n", time_now().c_str(), i, neurons.size());
      }
      neuron *src = neurons[i];

      // electrode contributions and (alternative) minimum distance fix
      src->electrode_contributions.reserve(electrodes.size());
      for (size_t k = 0; k < electrodes.size(); k++) {
        electrode *e = electrodes[k];
        double rik_squ = get_dist_squ(src, e);
        // this should not happen any more since we placed electrodes first.
        // if neuron too close to electrode, move neuron a bit
        if (rik_squ < r_min*r_min) {
          printf("neuron too close to electrode.\n");
          exit(-1);
          double dx = src->x - e->x;
          double dy = src->y - e->y;
          double a = dy/dx;
          double r = 1.01*pow(r_min*r_min/rik_squ, .5);
          double xn = e->x + dx*r;
          double yn = e->y + a*dx*r;
          src->x = xn-std::floor(xn/sys_size)*sys_size;
          src->y = yn-std::floor(yn/sys_size)*sys_size;
          rik_squ = get_dist_squ(src, e);
        }

        // coarse sampling contribution
        double irik = pow(rik_squ, -1./2.);
        src->electrode_contributions.push_back(irik);

        // subsampling closest neurons
        if (rik_squ < e->closest_neuron_distance_squ) {
          e->closest_neuron_distance_squ = rik_squ;
          e->closest_neuron = src;
        }
      }

      // activation probabilities based on distance
      size_t n_cout = src->outgoing.size();
      src->outgoing_probability.reserve(n_cout);
      double norm = 0.;
      for (size_t j = 0; j < n_cout; j++) {
        double dij_squ = get_dist_squ(src, src->outgoing[j]);
        double pij = exp(-dij_squ/sigma_squ);
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

  // neur_dist
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
        #ifdef NOCC
        // simple recurrent activations. if target already active, skip.
        if (!src->outgoing[j]->active
            && udis(rng) < src->outgoing_probability[j]) {
          activate_neuron(src->outgoing[j]);
        }
        #else
        // coalesence compensation
        // if random number says to activate, try current target.
        // if current target already active, try the next closest (farther)
        // connected neuron until successful or no more candidates
        if (udis(rng) < src->outgoing_probability[j]) {
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
        // only convert to lower precision float here, when storing,
        // not for individual neuron contributions to electrode signal
        real_t r = real_t(e->cs_signal);
        cs_histories[i].push_back(r);
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
    hdf5_write_scalar(h5file, "/meta/neur_dist",
                      sys->neur_dist, H5T_NATIVE_DOUBLE);
    hdf5_write_scalar(h5file, "/meta/elec_dist",
                      sys->elec_dist, H5T_NATIVE_DOUBLE);
    hdf5_write_scalar(h5file, "/meta/neuron_density",
                      sys->neuron_density, H5T_NATIVE_DOUBLE);
    hdf5_write_scalar(h5file, "/meta/d_max",
                      sys->d_max, H5T_NATIVE_DOUBLE);
    hdf5_write_scalar(h5file, "/meta/sigma",
                      sys->sigma, H5T_NATIVE_DOUBLE);
    hdf5_write_scalar(h5file, "/meta/m_micro",
                      sys->m_micro, H5T_NATIVE_DOUBLE);
    hdf5_write_scalar(h5file, "/meta/h_prob",
                      sys->h_prob, H5T_NATIVE_DOUBLE);

    // create data sets
    h5dset_cs = hdf5_create_appendable_nd(
                    h5file, "/data/coarse", H5T_CUSTOM_REAL, sys->num_elec,
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
      hdf5_append_nd(h5dset_cs, sys->cs_histories[i], H5T_CUSTOM_REAL, i,
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


// helper routine to iteratively find the desired activity by tuning drive
#ifdef FINDPAR
void find_parameters(size_t num_neur, double neur_dist,
                     size_t num_elec, double elec_dist,
                     size_t num_outgoing,
                     double m_micro, double sigma, double h_prob,
                     double delta_t,
                     double tau_target, double tau_prec,
                     double act_target, double act_prec,
                     size_t thrm_steps, size_t time_steps,
                     double dh, double dm
                    ) {

  double m_now = m_micro;
  double h_now = h_prob;
  bool found = false;


  while (!found) {
    printf("m = %e h = %e ", m_now, h_now);

    // silence stdout
    int old_stdout = dup(STDOUT_FILENO);
    FILE *nullout = fopen("/dev/null", "w");
    dup2(fileno(nullout), STDOUT_FILENO);

    asystem *sys = new asystem(num_neur, neur_dist, num_elec, elec_dist,
                               num_outgoing,
                               m_now, sigma, h_now,
                               time_steps);
    // restore stdout
    fclose(nullout);
    dup2(old_stdout, STDOUT_FILENO);
    close(old_stdout);
    setbuf(stdout, NULL);
    printf("created\n");

    for (size_t i = 0; i < thrm_steps; i++) {
      if(is_percent(i, thrm_steps, 1.))
        printf("therm: %02.1f%%\r", is_percent(i, thrm_steps, 1.));
      sys->update_step();
      sys->measure_step(i, false);
    }
    for (size_t i = 0; i < time_steps; i++) {
      if(is_percent(i, time_steps, 1.))
        printf("simul: %02.1f%%\r", is_percent(i, time_steps, 1.));
      sys->update_step();
      sys->measure_step(i, true);
    }

    double act, var, mlr, tau, act_hz;
    var = variance(sys->at_history, act);
    mlr = m_from_lin_regr(sys->at_hist_2d);

    tau = -delta_t/log(mlr);
    act_hz = act/sys->num_neur/delta_t*1000.;
    printf("last %.1e steps: act ~%4.2fHz, cv ~%.2f, ",
           double(sys->at_history.size()),
           act/sys->num_neur/delta_t*1000., sqrt(var)/act);
    printf("mlr ~%.5f, tau ~%.1fms\n\n",
           mlr, tau);

    found = true;
    if (fabs(tau - tau_target) > tau_prec && dm != 0) {
      if      (tau < tau_target) m_now += dm;
      else if (tau > tau_target) m_now -= dm;
      found = false;
    }
    if (fabs(act_hz - act_target) > act_prec && dh != 0) {
      if      (act_hz < act_target) h_now += dh;
      else if (act_hz > act_target) h_now -= dh;
      found = false;
    }
    delete sys;

  }
}
#endif


// ------------------------------------------------------------------ //
// main
// ------------------------------------------------------------------ //
int main(int argc, char *argv[]) {

  double time_steps   = 1e3;        // number of time steps
  double thrm_steps   = 1e3;        // thermalization steps before measuring
  size_t num_neur     = 256000;     // number of neurons
  size_t num_outgoing = 1000;       // average outgoing connections per neuron
  size_t num_elec     = 64;         // total number of electrodes
  double neur_dist    = 1.;         // inter-neuron (nearest-neigbour) distance
  double elec_dist    = 8.;         // electrode dist. [unit=nearestneur-dist]
  size_t seed         = 314;        // seed for the random number generator
  double m_micro      = .98;        // branching parameter applied locally
  double sigma        = 12.;        // eff. conn-length [unit=nearestneur-dist]
  double h            = 4e-5;       // probability for spontaneous activation
  double delta_t      = 2.;         // time step size [ms]
  double cache_size   = 1e5;        // [num time steps] before hist is written
  std::string path    = "";         // output path for results

  #ifdef FINDPAR
  double dh = 1e-5;                 // additive change to h
  double dm = 0e-5;                 // additive change to m
  double A  = 1.;                   // target activity in Hz
  double ta = .1;                   // activity tolerance in Hz
  #endif


  for (size_t i=0; i < argc; i++) {
    if (i+1 != argc) {
      if(std::string(argv[i])=="-T" ) time_steps     = atof(argv[i+1]);
      if(std::string(argv[i])=="-t" ) thrm_steps     = atof(argv[i+1]);
      if(std::string(argv[i])=="-N" ) num_neur       = atof(argv[i+1]);
      if(std::string(argv[i])=="-k" ) num_outgoing   = atof(argv[i+1]);
      if(std::string(argv[i])=="-e" ) num_elec       = atof(argv[i+1]);
      if(std::string(argv[i])=="-dn") neur_dist      = atof(argv[i+1]);
      if(std::string(argv[i])=="-de") elec_dist      = atof(argv[i+1]);
      if(std::string(argv[i])=="-s" ) seed           = atof(argv[i+1]);
      if(std::string(argv[i])=="-m" ) m_micro        = atof(argv[i+1]);
      if(std::string(argv[i])=="-g" ) sigma          = atof(argv[i+1]);
      if(std::string(argv[i])=="-h" ) h              = atof(argv[i+1]);
      if(std::string(argv[i])=="-c" ) cache_size     = atof(argv[i+1]);
      if(std::string(argv[i])=="-o" ) path           =      argv[i+1] ;

      #ifdef FINDPAR
      if(std::string(argv[i])=="-dh") dh             = atof(argv[i+1]);
      if(std::string(argv[i])=="-dm") dm             = atof(argv[i+1]);
      if(std::string(argv[i])=="-A")  A              = atof(argv[i+1]);
      if(std::string(argv[i])=="-ta") ta             = atof(argv[i+1]);
      #endif
    }
  }

  setbuf(stdout, NULL); // print direct without calling flush, also on cluster
  rng.seed(1000+seed);
  rng.discard(70000);

  #ifdef FINDPAR
  find_parameters(num_neur, neur_dist, num_elec, elec_dist,
                  num_outgoing,
                  m_micro, sigma, h,
                  delta_t,
                  400., 25.,
                  A, ta,
                  thrm_steps, time_steps,
                  dh, dm);
  return 0;
  #else
  if (path == "") {
    printf("specify output path with '-o'\n");
    return -1;
  }
  #endif


  asystem *sys = new asystem(num_neur, neur_dist, num_elec, elec_dist,
                             num_outgoing,
                             m_micro, sigma, h,
                             cache_size);

  exporter *out = new exporter(path, sys, seed);


  // hack to save some time by initializing to target cativity of 1Hz
  for (size_t i = 0; i < size_t(1.*sys->num_neur*delta_t/1000.); i++) {
    sys->activate_neuron(sys->neurons[size_t(udis(rng)*sys->num_neur)]);
  }

  printf("thermalizing for %.0e steps\n", thrm_steps);
  for (size_t i = 0; i < thrm_steps; i++) {
    sys->update_step();
    sys->measure_step(i, false);
  }


  printf("simulating for %.0e time steps\n", time_steps);
  for (size_t i = 0; i < size_t(time_steps); i++) {
    if(is_percent(i, size_t(time_steps), 5.) || have_passed_hours(6.)) {
      double act, var, mlr;
      var = variance(sys->at_history, act);
      mlr = m_from_lin_regr(sys->at_hist_2d);
      printf("\t%s, ~%2.0f%%, last %.1e steps: act ~%4.2fHz, cv ~%.2f, ",
             time_now().c_str(), is_percent(i, size_t(time_steps), 5.),
             double(sys->at_history.size()),
             act/sys->num_neur/delta_t*1000., sqrt(var)/act);
      printf("mlr ~%.5f, tau ~%.1fms\n",
             mlr, -delta_t/log(mlr));
      // out->write_histories();
    }
    sys->update_step();
    sys->measure_step(i, true);
    if(i%sys->cache_size==0) {
      out->write_histories();
      reset_act_hist(sys->at_hist_2d);
    }
  }
  out->write_histories();
  printf("done\n\n");
  out->finalize();

  return 0;
}
