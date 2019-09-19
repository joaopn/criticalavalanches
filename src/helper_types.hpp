#ifndef FPS_NEURO_TYPES
#define FPS_NEURO_TYPES

#include <vector>

class electrode;
class neuron {
 public:
  bool   active;
  double x;
  double y;
  double R_d;       // [um] radius of dendritic tree
  double l;         // [um] total length of the axon
  size_t id;        // needs to match index of neuron in the vector!
  std::vector< neuron * >    outgoing;
  std::vector< neuron * >    rejected;
  std::vector< double >      outgoing_probability;
  std::vector< electrode * > electrodes;
  std::vector< double >      electrode_contributions;

  neuron(size_t id_ = 0,
         double x_ = std::numeric_limits<double>::max(),
         double y_ = std::numeric_limits<double>::max()) {
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

#endif
