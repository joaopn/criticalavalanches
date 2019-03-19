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
#include <deque>
#include <random>
#include "assert.h"
#include <memory>         // std::unique_ptr
#include <zlib.h>         // compression

#define NDEBUG            // uncomment to enable assertion checks

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
  std::vector< neuron* >    outgoing;
  std::vector< electrode* > electrodes;
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
  size_t num_boxes;
  double sys_size;
  double elec_frac;
  double box_size;
  double delta_l;
  double delta_l_nn;
  double neuron_density;
  double outgoing_distance;
  double m_micro;
  double h_prob;
  std::vector< neuron* >    neurons;
  std::vector< electrode* > electrodes;

  // keep track of active neurons
  std::deque < neuron* > active_neurons;
  size_t num_active_new;
  size_t num_active_old;

  // signal recording
  std::vector< std::vector< double > > electrode_cs_histories;
  std::vector< std::vector< bool > >   electrode_ss_histories;

  // construct system
  asystem(double sys_size_ = 1., double elec_frac_ = .25,
    size_t num_neur_ = 144000, size_t num_outgoing_ = 1000,
    size_t num_elec_ = 100, double m_micro_ = 1.0, double h_prob_ = 0.1) {
    num_neur       = num_neur_;
    num_outgoing   = num_outgoing_;
    num_elec       = num_elec_;
    sys_size       = sys_size_;
    elec_frac      = elec_frac_;
    m_micro        = m_micro_;
    h_prob         = h_prob_;
    neuron_density = double(num_neur)/sys_size/sys_size;
    num_active_old = 0;
    num_active_new = 0;

    // analytic solution for average inter neuron distance delta_l
    delta_l = pow(sys_size, 3.)/6. * ( sqrt(2.) + log(1. + sqrt(2.)) );

    // analytic solution for average nearest-neighbour distance
    delta_l_nn = -exp(-neuron_density * M_PI)
      + erf(sqrt(neuron_density * M_PI))/2./sqrt(neuron_density);

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

    printf("system created\n");
    printf("\tnumber of neurons: %lu\n", num_neur);
    printf("\toutgoing connections per neuron: ~%lu\n", num_outgoing);
    // printf("\t\t(measured: %.2f)\n", mc);
    printf("\taverage distance between neurons: %.2e\n", delta_l);
    // printf("\t\t(measured: %.2e)\n", measure_avg_distance());
    printf("\taverage distance between nearest neighbours: %.2e\n",
      delta_l_nn);
    // printf("\t\t(measured: %.2e)\n", measure_avg_nearest_neighbour_distance());

    // create electrodes to spread over a frac of the system (each direction)
    size_t ne = 0;
    double de = sys_size*elec_frac/sqrt(double(num_elec));
    for (size_t i = 0; i < sqrt(num_elec); i++) {
      for (size_t j = 0; j < sqrt(num_elec); j++) {
        electrode *e = new electrode(j*de, i*de, ne);
        electrodes.push_back(e);
        ne += 1;
        // init electrode histories
        electrode_ss_histories.push_back(std::vector< bool >());
        electrode_ss_histories.back().reserve(1000);
        electrode_cs_histories.push_back(std::vector< double >());
        electrode_cs_histories.back().reserve(1000);
      }
    }

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
      printf("connecting: %lu/%lu\r", i, neurons.size());
    }
    printf("\33[2Kdone\n");
    return avg_connection_count/double(neurons.size());
  }

  void set_contributions_and_probabilities(double elec_distance) {
    printf("calculating contributions to each electrode\n");
    // calculate contributions to electrodes and probabilities to activate
    double w_squ = 2.*pow(4., 2.);
    // double m_micro = 1.1;
    // todo: this should NOT go with delta_l_nn but be fixed
    // double dik_min_squ = pow(.2*delta_l_nn, 2.);
    double dik_min_squ = pow(elec_distance/10., 2.);
    // if we do not ensure this, neurons will always be too close to some elec
    assert(dik_min_squ < pow(.5*elec_distance, 2.));
    for (size_t i = 0; i < neurons.size(); i++) {
      printf("contributions: %lu/%lu\r", i, neurons.size());
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
      // activation probabilities
      size_t n_cout = src->outgoing.size();
      src->outgoing_probability.reserve(n_cout);
      double norm = 0.;
      for (size_t j = 0; j < n_cout; j++) {
        double dij_squ = get_dist_squ(src, src->outgoing[j]);
        double pij = m_micro*exp(-dij_squ/w_squ);
        src->outgoing_probability.push_back(pij);
        norm += pij;
      }
      // normalize probabilities
      for (size_t j = 0; j < n_cout; j++) {
        src->outgoing_probability[j] /= norm;
      }
    }
    printf("\33[2Kdone\n");
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

  void update_step() {

    // double h_prob = 0.1;
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
      for (size_t j = 0; j < src->outgoing.size(); j++) {
        if (udis(rng) < src->outgoing_probability[j]
          && !src->outgoing[j]->active) {
          activate_neuron(src->outgoing[j]);
        }
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

    num_active_old = num_active_new;
    num_active_new = 0;
  }

  void measure_step() {
    assert(electrodes.size() == electrode_cs_histories.size());
    for (size_t i = 0; i < electrodes.size(); i++) {
      electrode *e = electrodes[i];
      electrode_cs_histories[i].push_back(e->cs_signal);
      electrode_ss_histories[i].push_back(e->ss_signal);
      e->cs_signal = 0.;
      e->ss_signal = false;
    }
  }
};


class exporter {
 public:
  std::vector< gzFile > gz_files;
  asystem *sys;

  exporter(std::string filepath, asystem *sys_, size_t seed) {
    sys = sys_;
    char filechar[2048];
    sprintf(filechar, "%s/N%06lu_m%.3f_s%05lu/",
      filepath.c_str(), sys->num_neur, sys->m_micro, seed);

    printf("exporting files to %s\n", filechar);
    fflush(stdout);
    system(("mkdir -p " + std::string(filechar)).c_str());

    for (size_t i = 0; i < sys->electrodes.size(); i++) {
      char temp[2048];
      sprintf(temp, "%s/elec%03lu.gz", filechar, i);
      gzFile zfile = gzopen(temp, "w");
      std::stringstream header;
      header << "#electrude_id=" << sys->electrodes[i]->id << std::endl;
      header << "#x=" << sys->electrodes[i]->x << std::endl;
      header << "#y=" << sys->electrodes[i]->y << std::endl;
      header << "#coarse_signal\t" << "ss_signal\n";
      gzprintf(zfile, header.str().c_str());
      gz_files.push_back(zfile);
    }
  }

  void finalize() {
    for (size_t i = 0; i < gz_files.size(); i++) {
      gzclose(gz_files[i]);
    }
  }

  // write cached electrode histores, assumes nested vectors of bool or double
  template <typename T1, typename T2>
  void write_histories(T1 &cs_histories, T2 &ss_histories) {
    assert(cs_histories.size() == gz_files.size());
    for (size_t i = 0; i < cs_histories.size(); i++) {
      for (size_t j = 0; j < cs_histories[i].size(); j++) {
        gzprintf(gz_files[i], "%e\t%lu\n",
          cs_histories[i][j], size_t(ss_histories[i][j]));
      }
      cs_histories[i].resize(0);
      ss_histories[i].resize(0);
    }
  }
  void write_histories() {
    write_histories(sys->electrode_cs_histories, sys->electrode_ss_histories);
  }

  // write system configuration. vectors with ptrs to electordes or neurons
  template <typename T1>
  void write_config(std::string filename,
    std::vector<T1> &config_to_save) {
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

int main(int argc, char* argv[]) {

  double sys_size     = 1.;         // sytem size
  double time_steps   = 1e3;        // number of time steps
  size_t num_neur     = 144000;     // number of neurons
  size_t num_outgoing = 1000;       // average outgoing connections per neuron
  size_t num_elec     = 100;        // total number of electrodes
  double elec_frac    = .25;        // electrodes cover part of system per dim
  size_t seed         = 314;        // seed for the random number generator
  double m_micro      = 1.0;        // branching parameter applied locally
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
      if(std::string(argv[i])=="-h") h              = atof(argv[i+1]);
      if(std::string(argv[i])=="-o") path           =      argv[i+1] ;
    }
  }

  if (path == "") {
    printf("specify output path with '-o'\n");
    return -1;
  }

  rng.seed(1000+seed);
  rng.discard(70000);

  asystem *sys = new asystem(sys_size, elec_frac, num_neur, num_outgoing,
    num_elec, m_micro, h);
  exporter *exp = new exporter(path, sys, seed);

  for (size_t i = 0; i < size_t(time_steps); i++) {
    printf("step %05lu, activity ~ %.3f\r",
      i, sys->num_active_old/double(num_neur));
    fflush(stdout);
    sys->update_step();
    sys->measure_step();
    if(i%100==0) exp->write_histories();
    if(10*i%size_t(time_steps)==0) printf("\33[2K~%.1f%%\n", i/time_steps*100);
  }
  printf("\33[2Kdone\n");
  exp->finalize();

  return 0;


  // ------------------------------------------------------------------ //
  // topology test
  // ------------------------------------------------------------------ //
  exp->write_config("/Users/paul/Desktop/neurons_2.dat", sys->neurons);
  exp->write_config("/Users/paul/Desktop/electrodes_2.dat", sys->electrodes);

  std::vector< neuron* > test(100, nullptr);
  for (size_t i = 0; i < sys->electrodes.size(); i++) {
    electrode *e = sys->electrodes[i];
    printf("%lu %e %e\n", e->id, e->x, e->y);
    neuron *n = e->closest_neuron;
    printf("\t%lu %e %e\n", n->id, n->x, n->y);
    test[i] = n;
  }
  exp->write_config("/Users/paul/Desktop/closest_2.dat", test);

  neuron *t = sys->neurons[size_t(udis(rng)*num_neur)];
  std::vector< neuron* > foo = {t};

  double d = 0.;
  for (size_t i = 0; i < t->outgoing.size(); i++) {
    d += t->outgoing_probability[i];
    printf("prob: %e %e\n", t->outgoing_probability[i], d);
  }

  for (size_t i = 0; i < t->electrode_contributions.size(); i++) {
    printf("contrib: %e\n", t->electrode_contributions[i]);
  }

  exp->write_config("/Users/paul/Desktop/conn_2.dat", foo);
  exp->write_config("/Users/paul/Desktop/conns_2.dat", t->outgoing);

  delete sys;

}
