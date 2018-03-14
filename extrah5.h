#ifndef __EXTRAH5_H__
#define __EXTRAH5_H__

#include <H5Cpp.h>
#include <vector>
#include <string>
#include <type_traits>

// My lovely collection of not very well written HDF5 utility functions
// Don't expect too much. Really.
// TODO: Do this properly. Would require probably require refactoring other stuff though.
namespace X5 {


  template <typename FR, typename TO>
  herr_t copy(FR &from, TO &to, const char *name) {
    return H5Ocopy(from.getId(), name, to.getId(), name, 0, 0);
  }


  template <typename T>
  htri_t exists(T location, const char *name) {
    // return H5Oexists_by_name(location.getId(), name, H5P_DEFAULT);
    return H5Lexists(location.getId(), name, H5P_DEFAULT);
  }


  template <typename FR, typename TO>
  void rcopy(FR &from, TO &to) {
    // No error reporting for you!
    int n_gp = from.getNumObjs();
    for (int i = 0; i < n_gp; ++i) {
      std::string name = from.getObjnameByIdx(i);
      copy(from, to, name.c_str());
    }
  }



  // Read a 1d dataset into a vector. typename T has to be default constructable
  template <typename T>
  std::vector<T> read_vector(H5::DataSet &ds, H5::PredType pt) {
    H5::DataSpace sp = ds.getSpace();
    hsize_t dims;
    sp.getSimpleExtentDims(&dims);
    std::vector<T> v(dims, T());
    hsize_t start = 0;
    sp.selectHyperslab(H5S_SELECT_SET, &dims, &start);
    H5::DataSpace sp_mem(1, &dims);
    ds.read(v.data(), pt, sp_mem, sp);
    return v;
  }
  // Special case of vector string
  template <typename T>
  std::vector<T> read_vector(H5::DataSet &ds, hsize_t strlen) {
    H5::StrType stringtype(0, strlen);
    H5::DataSpace sp = ds.getSpace();
    hsize_t dims;
    sp.getSimpleExtentDims(&dims);
    std::vector<T> v;
    for (size_t j = 0; j < dims; ++j) {
      std::string s;
      hsize_t start = j;
      hsize_t count = 1;
      sp.selectHyperslab(H5S_SELECT_SET, &count, &start);
      H5::DataSpace sp_mem(1, &count);
      ds.read(s, stringtype, sp_mem, sp);
      v.push_back(s);
    }
    return v;
  }


    // Read a 1d column in direction dim of an NDIMS dimensional dataset into a vector
  template <typename T, int NDIMS>
  std::vector<T> read_column(H5::DataSet &ds, hsize_t *start, hsize_t count, hsize_t dim, H5::PredType pt) {
    H5::DataSpace sp = ds.getSpace();
    hsize_t Ncount[NDIMS];
    for (int i = 0; i < NDIMS; ++i) {
      Ncount[i] = 1;
    }
    Ncount[dim] = count;
    std::vector<T> v(count, T());
    sp.selectHyperslab(H5S_SELECT_SET, Ncount, start);
    H5::DataSpace sp_mem(1, &count);
    ds.read(v.data(), pt, sp_mem, sp);
    return v;
  }
  // Special case where data is strings
  template <typename T, int NDIMS>
  std::vector<T> read_column(H5::DataSet &ds, hsize_t *start, hsize_t count, hsize_t dim, hsize_t strlen) {
    H5::StrType stringtype(0, strlen);
    H5::DataSpace sp = ds.getSpace();
    hsize_t Ncount[NDIMS];
    for (int i = 0; i < NDIMS; ++i) {
      Ncount[i] = 1;
    }
    std::vector<T> v;
    for (size_t i = 0; i < count; ++i) {
      start[dim] = i;
      sp.selectHyperslab(H5S_SELECT_SET, Ncount, start);
      std::string s;
      hsize_t c2[1] = {1};
      H5::DataSpace sp_mem(1, c2);
      ds.read(s, stringtype, sp_mem, sp);
      // if this is a char (is there a better way of testing?), grab just the first char
      if (std::is_fundamental<T>::value)
        v.push_back(s[0]);
      else
        v.push_back(T(*s.data()));
      // I have a strong suspicion this function is really awful. But it works, and I've got other things to work on.
      // Should probably be a template specialization instead
    }
    return v;
  }


  template <typename T, typename V, typename Y>
  void make_attr(T &location, const char *name, const V &value, Y type) {
      hsize_t dims[1] {1};
    H5::DataSpace sp(1, dims);
    H5::Attribute at = location.createAttribute(name, type, sp);
    at.write(type, &value);
  }


  template <typename FR, typename TO>
  void copy_attrs(FR &from, TO &to) {
    // PROBABLY WORKS FINE FOR INT AND DOUBLE
    // DEFINITELY NOT FINE FOR STRINGS ETC
    int n_attrs = from.getNumAttrs();
    for (int i = 0; i < n_attrs; ++i) {
      auto id = H5Aopen_idx(from.getId(), i);
      char name[20];
      H5Aget_name(id, 20, name);
      // std::cout << name << std::endl;
      auto a = from.openAttribute(name);
      auto dt = a.getDataType();
      double v; // Hopefully this is safe on all processors. (so long as both doubles and ints are 64 bit it should be fine)
      // FIXME Find a way of choosing the right type/allocating the right amount of memory
      a.read(dt, &v);
      make_attr(to, name, v, dt);
    }
  }


} // namespace X5 end

#endif
