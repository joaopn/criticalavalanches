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

// rng and distributions
std::mt19937 rng(316);
std::uniform_real_distribution<> udis(0, 1);
std::normal_distribution<> gdis(0, 1);
std::binomial_distribution<int> bdis(10, 0.5);

class electrode;
class neuron {
 public:
  bool active;
  double x;
  double y;
  size_t id;
  std::vector< neuron* > outgoing;
  std::vector< double > outgoing_probability;
  std::vector< electrode* > electrodes;
  std::vector< double > electrode_contributions;

  neuron(double x_, double y_, double id_ = 0) {
    active = false;
    x = x_;
    y = y_;
    id = id_;
  }
};

class electrode {
 public:
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
  std::vector< neuron* > neurons;
  std::vector< electrode* > electrodes;

  // construct system
  asystem(double sys_size_ = 1., double elec_frac_ = .25,
    size_t num_neur_ = 144000, size_t num_outgoing_ = 1000,
    size_t num_elec_ = 100) {
    num_neur = num_neur_;
    num_outgoing = num_outgoing_;
    num_elec = num_elec_;
    sys_size = sys_size_;
    elec_frac = elec_frac_;
    neuron_density = double(num_neur)/sys_size/sys_size;

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
    // connect_neurons_no_dc();

    printf("system created:\n");
    printf("\tnumber of neurons: %lu\n", num_neur);
    printf("\toutgoing connections per neuron: ~%lu\n", num_outgoing);
    printf("\t\t(measured: %.2f)\n", mc);
    printf("\taverage distance between neurons: %.2e\n", delta_l);
    // printf("\t\t(measured: %.2e)\n", measure_avg_distance());
    printf("\taverage distance between nearest neighbours: %.2e\n",
      delta_l_nn);
    printf("\t\t(measured: %.2e)\n", measure_avg_nearest_neighbour_distance());

    // create electrodes to spread over a quater of the system (each direction)
    size_t ne = 0;
    double de = sys_size*elec_frac/sqrt(double(num_elec));
    for (size_t i = 0; i < sqrt(num_elec); i++) {
      for (size_t j = 0; j < sqrt(num_elec); j++) {
        electrode *e = new electrode(j*de, i*de, ne);
        electrodes.push_back(e);
        ne += 1;
      }
    }

    // calculate contributions to electrodes and probabilities to activate
    double w_squ = 2.*pow(4., 2.);
    double m_loc = 1.01;
    // todo: this should NOT go with delta_l_nn but be fixed
    // double dik_min_squ = pow(.2*delta_l_nn, 2.);
    double dik_min_squ = pow(de/10., 2.);
    // if we do not ensure this, neurons will always be too close to some elec
    assert(dik_min_squ < pow(.5*de, 2.));
    for (size_t i = 0; i < neurons.size(); i++) {
      neuron *src = neurons[i];
      // electrode contributions (and minimum distance fix)
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
        double pij = m_loc*exp(-dij_squ/w_squ);
        src->outgoing_probability.push_back(pij);
        norm += pij;
      }
      // normalize probabilities
      norm /= double(n_cout);
      for (size_t j = 0; j < n_cout; j++) {
        src->outgoing_probability[j] /= norm;
      }
    }


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

  // check again later, if wanted
  // https://stackoverflow.com/questions/42061788/c-get-k-smallest-elementsindices-from-vector-with-ties
  void connect_neurons_no_dc() {
    for (size_t i = 0; i < neurons.size(); i++) {
      neuron *src = neurons[i];
      double shortest_distance = std::numeric_limits<double>::max();
      double min_distance      = std::numeric_limits<double>::min();

      printf("neurons: %lu/%lu\r", i, neurons.size());
      while (src->outgoing.size() < num_outgoing) {
        neuron *past = nullptr;

        for (size_t j = 0; j < neurons.size(); j++) {
          neuron *tar = neurons[j];
          double d = get_dist_squ(src, tar);
          if (tar == src) continue;
          if (d <= min_distance) continue;
          if (d <= shortest_distance) {
            shortest_distance = d;
            past = tar;
          }
        }

        src->outgoing.push_back(past);
        min_distance = shortest_distance;
        shortest_distance = std::numeric_limits<double>::max();
      }
    }
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

  template <typename T1>
  void save_config(std::string filename,
    std::vector<T1> &config_to_save) {
    std::ofstream file (filename, std::ofstream::out);
    file << "#L=" <<std::setprecision(0)<<std::fixed << sys_size <<std::endl;
    file << "#N=" << config_to_save.size() <<std::endl;
    file << "#x y" << std::endl;
    file << std::setprecision(7) << std::scientific;
    for (size_t i = 0; i<config_to_save.size(); i++) {
      file << config_to_save[i]->x << "\t" << config_to_save[i]->y;
      file << std::endl;
    }
    file.close();
    file.clear();
  }

  // overload for default argument of full system
  void save_config(std::string filename) {
    save_config(filename, neurons);
  }
};


int main(int argc, char* argv[]) {

  double sys_size     = 1.;
  double elec_frac    = .25;
  size_t num_neur     = 20000;
  size_t num_outgoing = 1000;
  size_t num_elec     = 100;

  asystem *sys = new asystem(sys_size, elec_frac, num_neur, num_outgoing, num_elec);
  sys->save_config("/Users/paul/Desktop/neurons_1.dat");
  sys->save_config("/Users/paul/Desktop/electrodes_1.dat", sys->electrodes);

  delete sys;
  // num_neur = 50000;
  // sys = new asystem(sys_size, elec_frac, num_neur, num_outgoing, num_elec);
  // sys->save_config("/Users/paul/Desktop/neurons_2.dat");
  // sys->save_config("/Users/paul/Desktop/electrodes_2.dat", sys->electrodes);

}
