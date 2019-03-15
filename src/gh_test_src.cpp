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
  double closest_neuron_distance;

  electrode(double x_, double y_, double id_ = 0) {
    x = x_;
    y = y_;
    id = id_;
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
  double box_size;
  double delta_l;
  double delta_l_nn;
  double neuron_density;
  double outgoing_distance;
  std::vector< neuron* > neurons;
  std::vector< electrode* > electrodes;

  std::vector< std::list<neuron *> > boxes;
  std::vector< std::vector<int> > neighbour_boxes;

  // construct system
  asystem(size_t num_neur_ = 144000, size_t num_outgoing_ = 1000,
    size_t num_elec_ = 100) {
    num_neur = num_neur_;
    num_outgoing = num_outgoing_;
    num_elec = num_elec_;
    sys_size = 1.;
    neuron_density = double(num_neur)/sys_size/sys_size;

    // analytic solution for average inter neuron distance delta_l
    delta_l = pow(sys_size, 3.)/6. * ( sqrt(2.) + log(1. + sqrt(2.)) );

    // analytic solution for average nearest-neighbour distance
    delta_l_nn = -exp(-neuron_density * M_PI)
      + erf(sqrt(neuron_density * M_PI))/2./sqrt(neuron_density);

    // choose r so that approx num_outgoing neighbours are accessible
    outgoing_distance = sqrt(sys_size*sys_size*num_outgoing/M_PI
        /double(num_neur));
    init_domain_decomp(outgoing_distance);

    // create neurons
    for (size_t i = 0; i < num_neur; i++) {
      neuron *n = new neuron(udis(rng)*sys_size, udis(rng)*sys_size, i);
      neurons.push_back(n);
    }

    // assign neurons to boxes of domain decomposition
    for (size_t i = 0; i < neurons.size(); i++) {
      set_box(neurons[i]);
      assert(get_box(neurons[i]) < num_boxes*num_boxes);
    }

    // create connections
    double mc = connect_neurons_using_interaction_radius(outgoing_distance);
    // connect_neurons_no_dc();
    // connect_neurons_using_connection_count();

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
    double de = sys_size/double(sqrt(num_elec));
    for (size_t i = 0; i < sqrt(num_elec); i++) {
      for (size_t j = 0; j < sqrt(num_elec); j++) {
        electrode *e = new electrode(j*de, i*de, ne);
        electrodes.push_back(e);
        ne += 1;
      }
    }

  }

  // ------------------------------------------------------------------ //
  // domain decomposition
  // ------------------------------------------------------------------ //
  inline int get_box(neuron *p) {
    int nx, ny;
    nx = std::floor(p->x/box_size);
    ny = std::floor(p->y/box_size);
    // printf("%lu %d %d %e %e\n", p->id, nx, ny, p->x, p->y);
    return nx + num_boxes*ny;
  }

  inline void set_box(neuron *p) {
    int my_box = get_box(p);
    boxes[my_box].push_front(p);
  }

  inline void remove_from_box(int my_box, neuron *p) {
    boxes[my_box].remove(p);
  }

  void init_domain_decomp(double target_box_size) {

    size_t nb = size_t(std::floor(sys_size/target_box_size));
    if (nb < 4) {
      printf("not enough boxes, skipping domain decomposition\n");
      num_boxes = 1;
      box_size = sys_size;
      boxes.push_back(std::list<neuron*> ());
      std::vector<int> boxes_to_check = {0};
      neighbour_boxes.push_back(boxes_to_check);
    } else {
      for (int i = 0; i<nb*nb; i++) {
        int ir = ((i%nb+1==nb)  ?     i+1-nb : i+1 );
        int il = ((i%nb-1<0)    ?     i-1+nb : i-1 );
        int iu = ((i-nb<0)      ? nb*nb+i-nb : i-nb);
        int id = ((i+nb>=nb*nb) ? i+nb-nb*nb : i+nb);
        int iur = ((iu%nb+1==nb) ? iu+1-nb : iu+1 );
        int iul = ((iu%nb-1<0)   ? iu-1+nb : iu-1 );
        int idr = ((id%nb+1==nb) ? id+1-nb : id+1 );
        int idl = ((id%nb-1<0)   ? id-1+nb : id-1 );
        std::vector<int> boxes_to_check = {i, ir, il, iu, id, iur, iul, idr, idl};
        neighbour_boxes.push_back(boxes_to_check);
        boxes.push_back(std::list<neuron*> ());
      }
      num_boxes = nb;
      box_size = sys_size/double(num_boxes);
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

  void connect_neurons_using_connection_count() {
    // connect to neighbours using domain decomp and connection count
    for (size_t i = 0; i < neurons.size(); i++) {
      printf("connecting: %lu/%lu\r", i, neurons.size());
      neuron *src = neurons[i];
      std::vector<int> boxes_to_check = neighbour_boxes[get_box(src)];
      double shortest_distance = std::numeric_limits<double>::max();
      double min_distance      = std::numeric_limits<double>::min();

      // printf("i: %lu\n", i);
      while (src->outgoing.size() < num_outgoing) {
        neuron *past = nullptr;
        for (auto const &b : boxes_to_check) {
          for (std::list<neuron*>::iterator tar = boxes[b].begin();
            tar != boxes[b].end(); ++tar) {
            // if (boxes[b].size() == 0) continue;
            if (src == (*tar)) continue;
            double d = get_dist_squ(src, *tar);
            if (d <= min_distance) continue;
            if (d <= shortest_distance) {
              shortest_distance = d;
              past = (*tar);
            }
          }
        }

        assert(past != nullptr && "num_outgoing < num_neur?");
        src->outgoing.push_back(past);
        min_distance = shortest_distance;
        shortest_distance = std::numeric_limits<double>::max();
      }
    }
    printf("\33[2Kdone\n");
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


  void init_electordes() {
  }
};


int main(int argc, char* argv[]) {

  asystem *sys = new asystem(100000);
  // sys->save_config("/Users/paul/Desktop/neurons.dat");
  // sys->save_config("/Users/paul/Desktop/electrodes.dat", sys->electrodes);

}
