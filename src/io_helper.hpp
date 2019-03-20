#include <time.h>         // clock_t, clock, CLOCKS_PER_SEC
#include <ctime>          // time_t
#include <cstdio>
#include <cstdarg>

#if __cplusplus < 201703L
  #include <memory>
#endif

// mimics printf into std::string
// https://codereview.stackexchange.com/questions/187183/create-a-c-string-using-printf-style-formatting
std::string format(const char *fmt, ...) {
  char buf[256];

  va_list args;
  va_start(args, fmt);
  const auto r = std::vsnprintf(buf, sizeof buf, fmt, args);
  va_end(args);

  if (r < 0)
    // conversion failed
    return {};

  const size_t len = r;
  if (len < sizeof buf)
    // we fit in the buffer
    return { buf, len };

  #if __cplusplus >= 201703L
  // C++17: Create a string and write to its underlying array
  std::string s(len, '\0');
  va_start(args, fmt);
  std::vsnprintf(s.data(), len+1, fmt, args);
  va_end(args);

  return s;
  #else
  // C++11 or C++14: We need to allocate scratch memory
  auto vbuf = std::unique_ptr<char[]>(new char[len+1]);
  va_start(args, fmt);
  std::vsnprintf(vbuf.get(), len+1, fmt, args);
  va_end(args);

  return { vbuf.get(), len };
  #endif
}


// returns the percentage i/of if it is modulo p, else 0. e.g. use in if()
inline double is_percent(size_t i, size_t of, double p) {
  if (p <= 0 || size_t(i*100./p)%of == 0) return i/double(of)*100.;
  else return 0.;
}


// check time since last calling following functions, e.g. true every 6 hours
clock_t past_hours = clock();
clock_t past_mins  = clock();
inline bool have_passed_hours(double hours) {
  clock_t now = clock();
  if (double(now - past_hours)/CLOCKS_PER_SEC/3600. >= hours) {
    past_hours = now;
    return true;
  }
  return false;
}

inline bool have_passed_mins(double mins) {
  clock_t now = clock();
  if (double(now - past_mins)/CLOCKS_PER_SEC/60. >= mins) {
    past_mins = now;
    return true;
  }
  return false;
}


// returns passed time since the provided clock_t
std::string time_since(const clock_t start) {
  clock_t end = clock();
  float seconds_used = double(end - start) / CLOCKS_PER_SEC;
  int days_used = int(seconds_used/86400.0);
  seconds_used = fmod(seconds_used,86400.0);
  int hours_used = int(seconds_used/3600.0);
  seconds_used = fmod(seconds_used,3600.0);
  int minutes_used = int(seconds_used/60.0);
  seconds_used = fmod(seconds_used,60.0);
  return format("%02dd %02dh %02dm %02ds",
      days_used, hours_used, minutes_used, int(seconds_used));
}

// returns current time as YYYY-MM_dd_HH-mm-ss
std::string time_now() {
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];
  time (&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer,sizeof(buffer),"%Y-%m-%d %H-%M-%S", timeinfo);
  return std::string(buffer);
}

