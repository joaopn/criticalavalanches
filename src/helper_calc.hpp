#ifndef FPS_CALC
#define FPS_CALC

#include <cstdio>
#include <cmath>
#include <vector>


// fill 2d activity histogram H[act_old][act_new] and resize if needed
template<typename T>
inline void update_act_hist(std::vector<T> &act_hist, size_t a_old, size_t a_new) {
  if (act_hist.size()        <= a_old) act_hist.resize(a_old+1);
  if (act_hist[a_old].size() <= a_new) act_hist[a_old].resize(a_new+1, 0);
  act_hist[a_old][a_new] += 1;
}

template<typename T>
void reset_act_hist(std::vector<T> &act_hist) {
  for (size_t i = 0; i < act_hist.size(); i++) {
    for (size_t j = 0; j < act_hist[i].size(); j++) {
      act_hist[i][j] = 0;
    }
  }
}

// takes a 2d vector of the form H[act_old][act_new]
// set i_start = 1 if NO inpust is assumed!
template<typename T>
double m_from_lin_regr(std::vector<T> &act_hist, size_t i_start = 0) {
  double a_new_av = 0;
  double a_old_av = 0;
  double count = 0;
  for (size_t i = i_start; i < act_hist.size(); i++) {
    for (size_t j = 0; j < act_hist[i].size(); j++) {
      a_old_av += i*act_hist[i][j];
      a_new_av += j*act_hist[i][j];
      count += act_hist[i][j];
    }
  }
  if (count != 0) a_new_av /= count;
  if (count != 0) a_old_av /= count;

  double cov = 0;
  double var = 0;
  double a_new;
  double a_old;
  for (size_t i = i_start; i < act_hist.size(); i++) {
    for (size_t j = 0; j < act_hist[i].size(); j++) {
      double x = act_hist[i][j];
      a_old = double(i);
      a_new = double(j);
      var += x*(a_old-a_old_av)*(a_old-a_old_av);
      cov += x*(a_new-a_new_av)*(a_old-a_old_av);
    }
  }

  return cov/var;
}

// calculate mean over vector
template<typename T>
double mean(std::vector< T > &vector) {
  double mean_value=0;
  for (size_t i = 0; i < vector.size(); i++) mean_value += vector[i];
  mean_value = mean_value/vector.size();
  return mean_value;
}

// variance (mean returned via reference in second argument, if provided)
template<typename T>
double variance(std::vector< T > &vector, double &mean_value) {
  mean_value = mean(vector);
  double variance = 0;
  for (unsigned int i=0; i<vector.size(); i++) {
    variance += (vector[i]-mean_value)*(vector[i]-mean_value);
  }
  variance = variance/static_cast<double>(vector.size()-1);
  return variance;
}

template<typename T>
double variance(std::vector< T > &vector) {
  double mean_value;
  variance(vector, mean_value);
}

#endif
