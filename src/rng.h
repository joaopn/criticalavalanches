// ------------------------------------------------------------------ //
// randum number generator and some distributions
// call as e.g. udis(rng)
// ------------------------------------------------------------------ //

#ifndef FPS_RNG
#define FPS_RNG

#include <random>
#include <cmath>

std::mt19937 rng(314);
std::uniform_real_distribution<> udis(0, 1);
std::normal_distribution<>       ndis(0, 1);
std::binomial_distribution<int>  bdis(10, 0.5);
std::weibull_distribution<>      rdis(2.0, 1.0);  // rayligh distribution

void init_rng(int seed=314) {
  rng.seed(seed);
  rng.discard(70000);
}

inline void set_udis_param(double min, double max) {
  udis.param(std::uniform_real_distribution<>::param_type(min, max));
}

inline void set_bdis_param(double a, double b) {
  bdis.param(std::binomial_distribution<>::param_type(a, b));
}

inline void set_ndis_param(double mean, double std) {
  ndis.param(std::normal_distribution<>::param_type(mean, std));
}

// convert parameters from weibull to rayligh distribution
inline void set_rdis_param(double std) {
  double a = 2.0;
  double b = sqrt(2.0)*std;
  rdis.param(std::weibull_distribution<>::param_type(a, b));
}

#endif
