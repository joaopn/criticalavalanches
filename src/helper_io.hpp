#ifndef FPS_IO
#define FPS_IO

#include <stdio.h>
#include <stdlib.h>
#include <sys/utsname.h>

#include <time.h>         // clock_t, clock, CLOCKS_PER_SEC
#include <ctime>          // time_t
#include <cstdio>
#include <cstdarg>
#include <hdf5.h>

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

// returns current time as YYYY-MM-dd HH-mm-ss
std::string time_now() {
  time_t rawtime;
  struct tm *timeinfo;
  char buffer[80];
  time (&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer,sizeof(buffer),"%Y-%m-%d %H-%M-%S", timeinfo);
  return std::string(buffer);
}


// ------------------------------------------------------------------ //
// hdf5 helper
// ------------------------------------------------------------------ //

herr_t status;

void hdf5_write_string(hid_t h5file, std::string name, std::string str) {
  const hsize_t dims[1] = {1};
  hid_t memtype, dspace, dset;

  memtype = H5Tcopy(H5T_C_S1);
  status  = H5Tset_size(memtype, str.size());

  dspace = H5Screate_simple(1, dims, nullptr);
  dset   = H5Dcreate(h5file, name.c_str(), memtype, dspace, H5P_DEFAULT,
                     H5P_DEFAULT,H5P_DEFAULT);
  status = H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, str.c_str());

  status = H5Dclose(dset);
  status = H5Sclose(dspace);
  status = H5Tclose(memtype);
}

void hdf5_write_platform_details(hid_t h5file) {
  struct utsname buffer;

  if (uname(&buffer) != 0) {
    printf("unable to retreive platform details\n");
    return;
  }

  hdf5_write_string(h5file, "/uname/system",  buffer.sysname);
  hdf5_write_string(h5file, "/uname/node",    buffer.nodename);
  hdf5_write_string(h5file, "/uname/release", buffer.release);
  hdf5_write_string(h5file, "/uname/version", buffer.version);
  hdf5_write_string(h5file, "/uname/machine", buffer.machine);

}

template<typename T>
void hdf5_write_scalar(hid_t h5file, std::string name, T &scalar,
                       hid_t h5dtype) {
  // data buffer with size 1
  const hsize_t dims[1] = {1};

  // create dataspace
  hid_t dspace = H5Screate_simple(1, dims, nullptr);
  hid_t dset   = H5Dcreate(h5file, name.c_str(), h5dtype, dspace, H5P_DEFAULT,
                           H5P_DEFAULT, H5P_DEFAULT);
  //think about converting if difference between T and h5dtype
  status = H5Dwrite(dset, h5dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &scalar);
  if (status < 0) {
    printf("Could not write successfully, aborting.");
    exit(-1);
  }
  status = H5Dclose(dset);
  status = H5Sclose(dspace);
}

// careful with casting: provide matching data types! vectors of bool dont work
template<typename T>
void hdf5_write_vector(hid_t h5file, std::string name, std::vector<T> &vector,
                       hid_t h5dtype) {
  // data buffer with length of input vector
  const hsize_t dims[1] = {vector.size()};
  T *buf = &vector[0];
  // create H5 dataspace
  hid_t dspace = H5Screate_simple(1, dims, nullptr);
  hid_t dset   = H5Dcreate(h5file, name.c_str(), h5dtype, dspace, H5P_DEFAULT,
                           H5P_DEFAULT, H5P_DEFAULT);
  //write array and close once done
  status = H5Dwrite(dset, h5dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  if (status < 0) {
    printf("Could not write successfully, aborting.");
    exit(-1);
  }
  status = H5Dclose(dset);
  status = H5Sclose(dspace);

}

// create 1D dataset for appending 1D vector
hid_t hdf5_create_appendable(hid_t h5file, std::string name, hid_t h5dtype,
                             size_t chunk_size = 1) {
  // init empty 1D array without size limit. chunks should match write size
  hsize_t       dims[1] = {0};
  hsize_t   max_dims[1] = {H5S_UNLIMITED};
  hsize_t chunk_dims[1] = {chunk_size};
  hid_t dspace = H5Screate_simple(1, dims, max_dims);
  hid_t dcpl   = H5Pcreate(H5P_DATASET_CREATE);
  // compress if available and enable chunking
  if (H5Zfilter_avail(H5Z_FILTER_DEFLATE)) status = H5Pset_deflate(dcpl, 3);
  status = H5Pset_chunk(dcpl, 1, chunk_dims);

  hid_t dset = H5Dcreate(h5file, name.c_str(), h5dtype, dspace, H5P_DEFAULT,
                         dcpl, H5P_DEFAULT);
  status = H5Sclose(dspace);
  status = H5Pclose(dcpl);

  return dset;
}

// append 1D array to 1D dataset
template<typename T>
void hdf5_append(hid_t dset, std::vector<T> &vector, hid_t h5dtype) {
  // determine dimension of added data, to set new extent
  hid_t     dspace = H5Dget_space(dset);
  hsize_t  dims[1] = {hsize_t(H5Sget_simple_extent_npoints(dspace))};
  hsize_t start[1] = {dims[0]};
  hsize_t strde[1] = {1};
  hsize_t count[1] = {vector.size()};
  hsize_t block[1] = {1};
  dims[0] = dims[0] + vector.size();
  // get adjusted space, where to store buf
  status = H5Dset_extent(dset, dims);
  dspace = H5Dget_space(dset);
  status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET,
                               start, strde, count, block);
  //write array and close once done
  T *buf = &vector[0];
  hsize_t out_dims[1] = {vector.size()};
  hid_t mspace = H5Screate_simple(1, out_dims, nullptr);
  status = H5Dwrite(dset, h5dtype, mspace, dspace, H5P_DEFAULT, buf);
  if (status < 0) {
    printf("Could not write successfully, aborting.");
    exit(-1);
  }
  status = H5Sclose(mspace);
  status = H5Sclose(dspace);
}

#define RANK 2
// create Nx1D dataset for appending 1D vectors.
hid_t hdf5_create_appendable_nd(hid_t h5file, std::string name, hid_t h5dtype,
                                size_t N = 1, size_t chunk_size = 1) {
  // int rank = 2;
  hsize_t       dims[RANK] = {N, 0};
  hsize_t   max_dims[RANK] = {H5S_UNLIMITED, H5S_UNLIMITED};
  hsize_t chunk_dims[RANK] = {1, chunk_size};
  hid_t dspace = H5Screate_simple(RANK, dims, max_dims);
  hid_t dcpl   = H5Pcreate(H5P_DATASET_CREATE);
  // compress if available
  if (H5Zfilter_avail(H5Z_FILTER_DEFLATE)) status = H5Pset_deflate(dcpl, 3);
  status = H5Pset_chunk(dcpl, RANK, chunk_dims);
  // status = H5Pset_fill_value(dcpl, h5dtype, &fill_val);
  hid_t dset = H5Dcreate(h5file, name.c_str(), h5dtype, dspace, H5P_DEFAULT,
                         dcpl, H5P_DEFAULT);
  status = H5Sclose(dspace);
  status = H5Pclose(dcpl);

  return dset;
}

// append 1D array to Nx1D dataset at N
template<typename T>
void hdf5_append_nd(hid_t dset, std::vector<T> &vector, hid_t h5dtype,
                    size_t N, long int offset = -1) {
  // fetch the current file space & determine dimension of current dataspace
  hid_t dspace = H5Dget_space(dset);
  hsize_t *dims = new hsize_t[RANK];
  hsize_t *max_dims = new hsize_t[RANK];
  status = H5Sget_simple_extent_dims(dspace, dims, max_dims);
  // printf("%d %d %d %d\n", dims[0], dims[1], max_dims[0], max_dims[1]);

  // workaround to write in each row without extending the whole dataset
  if (offset < 0) offset = dims[1];

  // determine dimension of added data, set new extent
  hsize_t start[RANK] = {N, size_t(offset)};
  hsize_t strde[RANK] = {1, 1};
  hsize_t count[RANK] = {1, vector.size()};
  hsize_t block[RANK] = {1, 1};
  if (N+1 > dims[0]) dims[0] = N+1;
  if (offset + vector.size() > dims[1]) dims[1] = offset + vector.size();
  status = H5Dset_extent(dset, dims);
  dspace = H5Dget_space(dset);
  status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET,
                               start, strde, count, block);
  // matching target space
  T *buf = &vector[0];
  hsize_t out_dims[2] = {1, vector.size()};
  hid_t mspace = H5Screate_simple(2, out_dims, nullptr);

  // write array
  status = H5Dwrite(dset, h5dtype, mspace, dspace, H5P_DEFAULT, buf);
  if (status < 0) {
    printf("Could not write successfully, aborting.");
    exit(-1);
  }
  status = H5Sclose(mspace);
  status = H5Sclose(dspace);
}

hid_t hdf5_create_file(std::string filepath) {
  if (filepath.rfind('/')!=std::string::npos) {
    std::string dirname = filepath;
    dirname.erase(dirname.rfind('/'));
    system(("mkdir -p " + dirname).c_str());
  }
  printf("exporting files to %s\n", filepath.c_str());

  // open file and create groups
  hid_t h5file = H5Fcreate(filepath.c_str(), H5F_ACC_TRUNC,
                           H5P_DEFAULT, H5P_DEFAULT);

  if(h5file < 0) exit(-1);

  H5Gcreate(h5file, "/data",       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gcreate(h5file, "/axons",      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gcreate(h5file, "/neurons",    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gcreate(h5file, "/electrodes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gcreate(h5file, "/meta",       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gcreate(h5file, "/uname",      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hdf5_write_platform_details(h5file);
  hdf5_write_string(h5file, "/uname/original_file_path",  filepath);

  return h5file;
}

// ------------------------------------------------------------------ //
// tests
// ------------------------------------------------------------------ //

void testnew() {
  hid_t h5file = H5Fcreate("/Users/paul/temp/ca_rep_1.hdf5", H5F_ACC_TRUNC,
                           H5P_DEFAULT, H5P_DEFAULT);


  H5Gcreate(h5file, "/general", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  double g = 9.5;
  hdf5_write_scalar(h5file, "/general/g", g, H5T_NATIVE_DOUBLE);

  std::vector<double> foo(1000, 0);
  for (size_t i = 0; i < 1000; i++) foo[i] = 10.*i;
  hdf5_write_vector(h5file, "/foo_vector", foo, H5T_NATIVE_DOUBLE);
  hdf5_write_vector(h5file, "/foo_vector_cast", foo, H5T_NATIVE_FLOAT);

  std::vector<unsigned> bar(10, 0);
  bar[3] = 1;
  hdf5_write_vector(h5file, "bar_vector", bar, H5T_NATIVE_UINT);

  hid_t dset_app = hdf5_create_appendable(h5file, "appended", H5T_NATIVE_DOUBLE);
  for (size_t i = 0; i < 4; i++) {
    std::vector<double> blub(100, i);
    hdf5_append(dset_app, blub, H5T_NATIVE_DOUBLE);
  }

  hid_t dset_nd = hdf5_create_appendable_nd(h5file, "nd_appended",
                                            H5T_NATIVE_DOUBLE, 4, size_t(10));

  size_t j = 0;
  for (size_t i = 0; i < 4; i++) {
    std::vector<double> blub(size_t(10), j);
    j += 1;
    // default offset -1 is if only one row is written
    hdf5_append_nd(dset_nd, blub, H5T_NATIVE_DOUBLE, 0);
  }

  for (size_t i = 0; i < 4; i++) {
    std::vector<double> blub(size_t(10), j);
    j += 1;
    hdf5_append_nd(dset_nd, blub, H5T_NATIVE_DOUBLE, 5, i*blub.size());
  }

  for (size_t i = 3; i < 5; i++) {
    std::vector<double> blub(size_t(12), j);
    j += 1;
    hdf5_append_nd(dset_nd, blub, H5T_NATIVE_DOUBLE, 2, i*blub.size());
  }
}

#endif
