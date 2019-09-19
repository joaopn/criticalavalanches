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
#include "helper_types.hpp"
#include "electrode_sampling.cpp"
#ifdef TPORLANDI
  #include "topology_orlandi.cpp"
#elif defined TPGAUSS
  #include "topology_local_gauss.cpp"
#elif defined TPRANDOM
  #include "topology_random.cpp"
#endif

int main(int argc, char *argv[]) {

  double time_steps   = 1e5;        // number of time steps
  double thrm_steps   = 1e3;        // thermalization steps before measuring
  size_t num_neur     = 160000;     // number of neurons
  size_t num_outgoing = 1000;       // average outgoing connections per neuron
  size_t num_elec     = 64;         // total number of electrodes
  double neur_dist    = 50.;        // [um] inter-neur. (nearest-neigbour) dist
  double elec_dist    = 8.;         // [nearestneur-dist] electrode dist.
  double gamma        = -1.0;       // decay exp of dist-dep. contrib. to elec.
  size_t seed         = 314;        // seed for the random number generator
  double m_micro      = .98;        // branching parameter applied locally
  double sigma        = 6.;         // [nearestneur-dist] eff. conn-length
  double h            = 4e-5;       // probability for spontaneous activation
  double delta_t      = 2.;         // [ms] time step size
  double cache_size   = 1e5;        // [num time steps] before writing history
  size_t write_detail = 0;          // include topol. details (conn. matrix)
  std::string path    = "";         // output path for results


  for (size_t i=0; i < argc; i++) {
    if (i+1 != argc) {
      if(std::string(argv[i])=="-T" )  time_steps    = atof(argv[i+1]);
      if(std::string(argv[i])=="-t" )  thrm_steps    = atof(argv[i+1]);
      if(std::string(argv[i])=="-N" )  num_neur      = atof(argv[i+1]);
      if(std::string(argv[i])=="-k" )  num_outgoing  = atof(argv[i+1]);
      if(std::string(argv[i])=="-e" )  num_elec      = atof(argv[i+1]);
      if(std::string(argv[i])=="-dn")  neur_dist     = atof(argv[i+1]);
      if(std::string(argv[i])=="-de")  elec_dist     = atof(argv[i+1]);
      if(std::string(argv[i])=="-s" )  seed          = atof(argv[i+1]);
      if(std::string(argv[i])=="-m" )  m_micro       = atof(argv[i+1]);
      if(std::string(argv[i])=="-g" )  sigma         = atof(argv[i+1]);
      if(std::string(argv[i])=="-ga")  gamma         = atof(argv[i+1]);
      if(std::string(argv[i])=="-h" )  h             = atof(argv[i+1]);
      if(std::string(argv[i])=="-c" )  cache_size    = atof(argv[i+1]);
      if(std::string(argv[i])=="-x" )  write_detail  = atof(argv[i+1]);
      if(std::string(argv[i])=="-o" )  path          =      argv[i+1] ;
    }
  }

  setbuf(stdout, NULL); // print direct without calling flush, also on cluster
  if (path == "") {
    printf("specify output path with '-o'\n");
    return -1;
  }
  if (write_detail == 1 && num_neur > 20000) {
    printf("writing detailed output ('-x 1') is discouraged for large ");
    printf("numbers of neurons > 20000\n\toverwrite with '-x 2'\n");
    return -1;
  }

  hid_t h5file;
  double sys_size, xy_offset;
  size_t cc;
  std::vector<neuron *> neurons;

  init_rng(1000+seed);
  h5file = hdf5_create_file(path);
  hdf5_write_scalar(h5file, "/meta/seed", seed, H5T_NATIVE_HSIZE);

  // place electrode array in the middle of the culture
  sys_size  = 2.*sqrt(double(num_neur))*neur_dist;
  xy_offset = .5*sys_size - .5*sqrt(num_elec)*elec_dist*neur_dist;
  auto sam  = electrode_sampling();
  sam.par.NE     = num_elec;
  sam.par.d_E    = neur_dist*elec_dist;  // [um]
  sam.par.d_zone = neur_dist/5.;         // [um]
  sam.par.gamma  = gamma;
  sam.par.cache  = cache_size;
  sam.init_electrodes(xy_offset);
  sam.init_writing_hdf5(h5file);
  sam.print_parameters();

  // create neuron topology and avoid electrode positions
  #ifdef TPGAUSS
  auto tpl = topology_local_gauss();
  tpl.set_N_and_d_N(num_neur, neur_dist);
  tpl.par.K   = num_outgoing;
  tpl.par.std = sigma;
  #elif defined TPORLANDI
  auto tpl = topology_orlandi();
  tpl.par.rho = 1000.*1000./4./neur_dist/neur_dist; // [1/mm2]
  tpl.par.L   = sys_size;                           // [um]
  tpl.par.N   = num_neur;
  tpl.par.d_N = neur_dist;
  #elif defined TPRANDOM
  auto tpl = topology_random();
  tpl.par.rho = 1000.*1000./4./neur_dist/neur_dist; // [1/mm2]
  tpl.par.L   = sys_size;                           // [um]
  tpl.par.N   = num_neur;
  tpl.par.K   = num_outgoing;
  tpl.par.std = sigma;
  #endif
  if (write_detail > 0)
    tpl.init(neurons, h5file, sam.electrodes, pow(sam.par.d_zone, 2.));
  else
    tpl.init(neurons, -1, sam.electrodes, pow(sam.par.d_zone, 2.));
  tpl.print_parameters();

  // check that electrode array is within the culture
  for (size_t i = 0; i < sam.electrodes.size(); i++) {
    electrode *e = sam.electrodes[i];
    if (e->x < 0 || e->x > tpl.par.L || e->y < 0 || e->y > tpl.par.L) {
      printf("electrode array exceeds culture dimension\n");
      return -1;
    }
  }

  // set contribution strength of spikes to electrode signals
  sam.set_contributions(neurons, tpl);

  // set up dynamics and adjust recurrent activation probability to m
  auto dyn = dynamic_branching(neurons, sam.electrodes, m_micro, h);
  dyn.scale_activation_probabilities(neurons);

  // write details to file
  sam.write_sampling_details(h5file);
  dyn.write_dynamic_details(h5file);
  tpl.write_topology_details(h5file);
  if (write_detail > 0) {
    tpl.write_connectivty_matrix(neurons, h5file);
    tpl.write_neuron_details(neurons, h5file);
  }


  // ------------------------------------------------------------------ //
  // simulation
  // ------------------------------------------------------------------ //

  set_udis_param(0.0, 1.0);   // uniform dist min=0 max=1

  // hack to save some time by initializing to target activity of 1Hz
  for (size_t i = 0; i < size_t(1.*neurons.size()*delta_t/1000.); i++) {
    dyn.activate_neuron(neurons.at(size_t(udis(rng)*neurons.size())));
  }

  printf("thermalizing for %.0e steps\n", thrm_steps);
  for (size_t i = 0; i < thrm_steps; i++) {
    dyn.update_step();
    sam.measure_step(i, false);
  }

  size_t num_active;
  printf("simulating for %.0e time steps\n", time_steps);
  for (size_t i = 0; i < size_t(time_steps); i++) {
    if(is_percent(i, size_t(time_steps), 5.) || have_passed_hours(6.)) {
      double act, var, mlr;
      var = variance(sam.at_history, act);
      mlr = m_from_lin_regr(sam.at_hist_2d);
      printf("\t%s, ~%2.0f%%, last %.1e steps: act ~%4.2fHz, cv ~%.2f, ",
             time_now().c_str(), is_percent(i, size_t(time_steps), 5.),
             double(sam.at_history.size()),
             act/double(neurons.size())/delta_t*1000., sqrt(var)/act);
      printf("mlr ~%.5f, tau ~%.1fms\n",
             mlr, -delta_t/log(mlr));
    }
    dyn.update_step(num_active);
    sam.measure_step(i, num_active, true);
    if(i%sam.par.cache==0) {
      sam.write_histories();
      reset_act_hist(sam.at_hist_2d);
    }
  }
  sam.write_histories();

  for (size_t i = 0; i < neurons.size(); i++) delete neurons[i];
  H5Fclose(h5file);

  return 0;
}
