/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/HDF5pp

================================================================================================= */

#ifndef HDF5PP_H
#define HDF5PP_H

// ==================================== PREPROCESSOR DIRECTIVES ====================================

// -------------------------------- load libraries (conditionally) ---------------------------------

// basic include
#include <fstream>
#include "H5Cpp.h"
#include <vector>
#include <assert.h>

// optionally enable plug-in Eigen and load the library
#ifdef EIGEN_WORLD_VERSION
#define HDF5PP_EIGEN
#endif

#ifdef HDF5PP_EIGEN
#include <Eigen/Eigen>
#endif

// optionally enable plug-in cppmat and load the library
#ifdef CPPMAT_WORLD_VERSION
#define HDF5PP_CPPMAT
#endif

#ifdef HDF5PP_CPPMAT
#include <cppmat/cppmat.h>
#endif

// optionally enable plug-in xtensor and load the library
#ifdef XTENSOR_VERSION_MAJOR
#define HDF5PP_XTENSOR
#endif

#ifdef HDF5PP_XTENSOR
#include <xtensor/xarray.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xeval.hpp>
#include <xtensor/xexpression.hpp>
#include <xtensor/xio.hpp>
#endif

// -------------------------------------- version information --------------------------------------

#define HDF5PP_WORLD_VERSION 0
#define HDF5PP_MAJOR_VERSION 1
#define HDF5PP_MINOR_VERSION 4

#define HDF5PP_VERSION_AT_LEAST(x,y,z) \
  (HDF5PP_WORLD_VERSION>x || (HDF5PP_WORLD_VERSION>=x && \
  (HDF5PP_MAJOR_VERSION>y || (HDF5PP_MAJOR_VERSION>=y && \
                              HDF5PP_MINOR_VERSION>=z))))

#define HDF5PP_VERSION(x,y,z) \
  (HDF5PP_WORLD_VERSION==x && \
   HDF5PP_MAJOR_VERSION==y && \
   HDF5PP_MINOR_VERSION==z)

// ---------------------------- contain everything in the namespace H5p ----------------------------

namespace H5p {

// ======================================= SUPPORT FUNCTIONS =======================================

template<typename T> inline H5::PredType getType();

// ================================== CLASS DEFINTION (OVERVIEW) ===================================

class File
{
private:
  H5::H5File  m_file;
  std::string m_fname;
  bool        m_autoflush;

public:

  // constructor
  // -----------

  File() = default;

  File(const std::string &fname, const std::string &mode="w", bool autoflush=true);

  // support functions
  // -----------------

  // return the filename
  std::string fname() const;

  // flush all buffers associated with a file to disk
  // NB if 'autoflush==true' you don't need to call this function, all 'write' functions call it
  void flush();

  // check if a path exists (is a group or a dataset)
  bool exists(const std::string &path) const;

  // create a group
  // NB there is usually no need to call this function, all 'write' functions call it
  void createGroup(std::string path);

  // unlink a path
  // WARNING the space in the file may not be freed, use: $ h5repack file1 file2
  void unlink(std::string path);

  // read the shape of the data
  std::vector<size_t> shape(std::string path);

  // read the shape of the data along a specific dimension
  size_t shape(std::string path, size_t i);

  // read the size of the data (total number of entries, '== prod(shape(path))')
  size_t size(std::string path);

  // (advanced) read the size of the data in an opened dataset or data-space
  size_t size(const H5::DataSet   &dataset  );
  size_t size(const H5::DataSpace &dataspace);

  // (advanced) read the shape of the data in an opened dataset or data-space
  std::vector<size_t> shape(const H5::DataSet   &dataset  );
  std::vector<size_t> shape(const H5::DataSpace &dataspace);

  // (advanced) check if an opened dataset has the exact precision of the template type
  template<typename T>
  bool correct_presision(const H5::DataSet &dataset);

  // read from file
  // --------------

  // read as specific data-type, e.g. double, Eigen<double,...>, cppmat::array<double>, ...
  // NB double, int, etc. -> data can contain only exactly one entry
  //    std::string       -> data can contain only a string
  template<typename T>
  T read(std::string path);

  // read array component as specific scalar data-type, e.g. double, int, ...
  template<typename T>
  T read(std::string path, size_t index);

  // (advanced) read scalar of arbitrary type from a dataset containing exactly one entry
  template<typename T>
  T read_scalar(std::string path, const H5::PredType& HT);

  // (advanced) read scalar of arbitrary type from a dataset (scalar, or of rank 1)
  template<typename T>
  T read(std::string path, const H5::PredType& HT, size_t index);

  // (advanced) read "std::vector" of arbitrary type from a dataset of arbitrary rank
  template<typename T>
  std::vector<T> read_vector(std::string path, const H5::PredType& HT);

  // write to file
  // -------------

  // write "std::string" to string dataset
  void write(std::string path, std::string data);

  // write scalar to scalar dataset (non-extendable)
  void write(std::string path, int    data);
  void write(std::string path, size_t data);
  void write(std::string path, float  data);
  void write(std::string path, double data);

  // write scalar as (part of) an extendable dataset of rank 1
  void write(std::string path, int    data, size_t index, int    fillval=0  , size_t chunk_size=64000);
  void write(std::string path, size_t data, size_t index, size_t fillval=0  , size_t chunk_size=64000);
  void write(std::string path, float  data, size_t index, float  fillval=0.0, size_t chunk_size=64000);
  void write(std::string path, double data, size_t index, double fillval=0.0, size_t chunk_size=64000);

  // write "std::vector" to a dataset of arbitrary shape
  void write(std::string path, const std::vector<int>    &data,const std::vector<size_t> &shape={});
  void write(std::string path, const std::vector<size_t> &data,const std::vector<size_t> &shape={});
  void write(std::string path, const std::vector<float>  &data,const std::vector<size_t> &shape={});
  void write(std::string path, const std::vector<double> &data,const std::vector<size_t> &shape={});

  // (advanced) write scalar of arbitrary type to a dataset containing exactly one entry
  template<typename T>
  void write(std::string path, T data, const H5::PredType& HT);

  // (advanced) write scalar of arbitrary type as (part of) an extendable dataset of rank 1
  template<typename T>
  void write(std::string path, T data, const H5::PredType& HT, size_t index, T fill_val,
    size_t chunk_size);

  // (advanced) write array or any type and of arbitrary shape or rank
  template<typename T>
  void write(std::string path, const T *input, const H5::PredType& HT,
    const std::vector<size_t> &shape);

  // (advanced) write std::vector of arbitrary type to a dataset of arbitrary rank
  template<typename T>
  void write(std::string path, const std::vector<T> &data, const H5::PredType& HT,
    const std::vector<size_t> &shape);

  // overwrite to file
  // -----------------

  // overwrite scalar to scalar dataset (non-extendable)
  void overwrite(std::string path, int    data);
  void overwrite(std::string path, size_t data);
  void overwrite(std::string path, float  data);
  void overwrite(std::string path, double data);

  // overwrite "std::vector" to a dataset of arbitrary shape
  void overwrite(std::string path, const std::vector<int>    &data,const std::vector<size_t> &shape={});
  void overwrite(std::string path, const std::vector<size_t> &data,const std::vector<size_t> &shape={});
  void overwrite(std::string path, const std::vector<float>  &data,const std::vector<size_t> &shape={});
  void overwrite(std::string path, const std::vector<double> &data,const std::vector<size_t> &shape={});

  // (advanced) overwrite scalar of arbitrary type to a dataset containing exactly one entry
  template<typename T>
  void overwrite(std::string path, T data, const H5::PredType& HT);

  // (advanced) overwrite array or any type and of arbitrary shape or rank
  template<typename T>
  void overwrite(std::string path, const T *input, const H5::PredType& HT,
    const std::vector<size_t> &shape);

  // (advanced) overwrite std::vector of arbitrary type to a dataset of arbitrary rank
  template<typename T>
  void overwrite(std::string path, const std::vector<T> &data, const H5::PredType& HT,
    const std::vector<size_t> &shape);

  // plugin: Eigen
  // -------------

  #ifdef HDF5PP_EIGEN

  // write column to dataset of rank 1
  void write(std::string path, const Eigen::Matrix<int   ,Eigen::Dynamic,1,Eigen::ColMajor> &data);
  void write(std::string path, const Eigen::Matrix<size_t,Eigen::Dynamic,1,Eigen::ColMajor> &data);
  void write(std::string path, const Eigen::Matrix<float ,Eigen::Dynamic,1,Eigen::ColMajor> &data);
  void write(std::string path, const Eigen::Matrix<double,Eigen::Dynamic,1,Eigen::ColMajor> &data);

  // write matrix to dataset of rank 2
  void write(std::string path,
    const Eigen::Matrix<int  ,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &data);

  void write(std::string path,
    const Eigen::Matrix<size_t,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &data);

  void write(std::string path,
    const Eigen::Matrix<float ,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &data);

  void write(std::string path,
    const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &data);

  // (advanced) write column of arbitrary type to dataset of rank 1
  template<typename T>
  void write(std::string path,
    const Eigen::Matrix<T,Eigen::Dynamic,1,Eigen::ColMajor> &data,
    const H5::PredType& HT);

  // (advanced) write matrix of arbitrary type to dataset of rank 2
  template<typename T>
  void write(std::string path,
    const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &data,
    const H5::PredType& HT);

  // overwrite column to dataset of rank 1
  void overwrite(std::string path, const Eigen::Matrix<int   ,Eigen::Dynamic,1,Eigen::ColMajor> &data);
  void overwrite(std::string path, const Eigen::Matrix<size_t,Eigen::Dynamic,1,Eigen::ColMajor> &data);
  void overwrite(std::string path, const Eigen::Matrix<float ,Eigen::Dynamic,1,Eigen::ColMajor> &data);
  void overwrite(std::string path, const Eigen::Matrix<double,Eigen::Dynamic,1,Eigen::ColMajor> &data);

  // overwrite matrix to dataset of rank 2
  void overwrite(std::string path,
    const Eigen::Matrix<int  ,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &data);

  void overwrite(std::string path,
    const Eigen::Matrix<size_t,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &data);

  void overwrite(std::string path,
    const Eigen::Matrix<float ,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &data);

  void overwrite(std::string path,
    const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &data);

  // (advanced) overwrite column of arbitrary type to dataset of rank 1
  template<typename T>
  void overwrite(std::string path,
    const Eigen::Matrix<T,Eigen::Dynamic,1,Eigen::ColMajor> &data,
    const H5::PredType& HT);

  // (advanced) overwrite matrix of arbitrary type to dataset of rank 2
  template<typename T>
  void overwrite(std::string path,
    const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &data,
    const H5::PredType& HT);

  // (advanced) read data of arbitrary type to Eigen column
  template<typename T>
  Eigen::Matrix<T,Eigen::Dynamic,1,Eigen::ColMajor> read_eigen_column(std::string path,
    const H5::PredType& HT);

  // (advanced) read data of arbitrary type to Eigen matrix
  template<typename T>
  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> read_eigen_matrix(std::string path,
    const H5::PredType& HT);

  #endif

  // plugin: cppmat
  // --------------

  #ifdef HDF5PP_CPPMAT

  // write nd-array to dataset of matching rank
  void write(std::string path, const cppmat::array<int   > &data);
  void write(std::string path, const cppmat::array<size_t> &data);
  void write(std::string path, const cppmat::array<float > &data);
  void write(std::string path, const cppmat::array<double> &data);

  // (advanced) write nd-array of arbitrary type to dataset of matching rank
  template<typename T>
  void write(std::string path, const cppmat::array<T> &data, const H5::PredType& HT);

  // overwrite nd-array to dataset of matching rank
  void overwrite(std::string path, const cppmat::array<int   > &data);
  void overwrite(std::string path, const cppmat::array<size_t> &data);
  void overwrite(std::string path, const cppmat::array<float > &data);
  void overwrite(std::string path, const cppmat::array<double> &data);

  // (advanced) overwrite nd-array of arbitrary type to dataset of matching rank
  template<typename T>
  void overwrite(std::string path, const cppmat::array<T> &data, const H5::PredType& HT);

  // (advanced) read data of arbitrary type to cppmat::array
  template<typename T>
  cppmat::array<T> read_cppmat_array(std::string path, const H5::PredType& HT);

  #endif

  // plugin: xtensor <TAKEN FROM https://github.com/QuantStack/xtensor-io>
  // the feature from this library will be moved there and to HighFive,
  // this is just a temporary fix
  // ---------------------------------------------------------------------

  #ifdef HDF5PP_XTENSOR

  // write nd-array to dataset of matching rank
  template<class E> void write(std::string path, const xt::xexpression<E> &data);

  // overwrite nd-array to dataset of matching rank
  template<class E> void overwrite(std::string path, const xt::xexpression<E> &data);

  // read "xarray" of arbitrary type
  template<class T> auto xread(const std::string& path);

  // read "xtensor" of arbitrary type
  template<class T, size_t N> auto xread(const std::string& path);

  // (advanced) generic read
  template<class T> T xread_impl(const std::string& path, const H5::PredType& HT);

  #endif
};

// ======================================= SUPPPORT FUNCTION =======================================

template<> inline H5::PredType getType<int   >() { return H5::PredType::NATIVE_INT;    }
template<> inline H5::PredType getType<size_t>() { return H5::PredType::NATIVE_HSIZE;  }
template<> inline H5::PredType getType<float >() { return H5::PredType::NATIVE_FLOAT;  }
template<> inline H5::PredType getType<double>() { return H5::PredType::NATIVE_DOUBLE; }

// ========================================= CONSTRUCTORS ==========================================

inline File::File(const std::string &name, const std::string &mode_, bool autoflush)
{
  // copy mode
  std::string mode = mode_;

  // copy filename
  m_fname = name;

  // check if the file exists
  if ( mode == "r" )
  {
     // - find file
    std::ifstream infile(m_fname);
    // - throw error if file does not exist
    if ( ! infile.good() ) std::runtime_error("HDF5pp: file does not exist ('"+name+"')");
  }

  // check if file exists, otherwise set write mode to "w"
  if ( mode == "a" or mode == "r+" )
  {
    // - find file
    std::ifstream infile(m_fname);
    // - change write mode if file does not exist
    if ( ! infile.good() ) mode = "w";
  }

  // open file
  if      ( mode == "r"         ) m_file = H5::H5File(m_fname.c_str(),H5F_ACC_RDONLY);
  else if ( mode == "w"         ) m_file = H5::H5File(m_fname.c_str(),H5F_ACC_TRUNC );
  else if ( mode == "a" or "r+" ) m_file = H5::H5File(m_fname.c_str(),H5F_ACC_RDWR  );
  else throw std::runtime_error("HDF5pp: unknown mode '"+mode+"'");

  // store flush settings
  m_autoflush = autoflush;
}

// ======================================= SUPPORT FUNCTIONS =======================================

// ---------------------------------------- return filename ----------------------------------------

inline std::string File::fname() const
{
  return m_fname;
}

// ------------------------------------------ flush file -------------------------------------------

inline void File::flush()
{
  m_file.flush(H5F_SCOPE_GLOBAL);
}

// -------------------------- check if path exists (is group or dataset) --------------------------

inline bool File::exists(const std::string &path) const
{
  // find first "/"
  size_t idx = path.find("/");

  // loop over all groups
  while ( true )
  {
    // - terminate if all "/" have been found
    if ( std::string::npos == idx ) break;
    // - create group if needed
    if ( idx > 0 )
    {
      // -- get group name
      std::string name(path.substr(0,idx));
      // -- create if needed
      if ( !m_file.exists(name.c_str()) ) return false;
    }
    // - proceed to next "/"
    idx = path.find("/",idx+1);
  }

  return m_file.exists(path.c_str());
}

// ---------------------------------------- create a group -----------------------------------------

inline void File::createGroup(std::string path)
{
  // find first "/"
  size_t idx = path.find("/");

  // loop over all groups
  while ( true )
  {
    // - terminate if all "/" have been found
    if ( std::string::npos == idx ) return;
    // - create group if needed
    if ( idx > 0 )
    {
      // -- get group name
      std::string name(path.substr(0,idx));
      // -- create if needed
      if ( !m_file.exists(name.c_str()) )
        H5::Group group = m_file.createGroup(name.c_str());
    }
    // - proceed to next "/"
    idx = path.find("/",idx+1);
  }
}

// ----------------------------------------- unlink a path -----------------------------------------

inline void File::unlink(std::string path)
{
  m_file.unlink(path.c_str());
}

// ------------------------ read size of the data (total number of entries) ------------------------

inline size_t File::size(std::string path)
{
  // check existence of path
  if ( ! exists(path) )
    throw std::runtime_error("HDF5pp::size: dataset not found ('"+path+"')");

  // return size
  return static_cast<size_t>(m_file.openDataSet(path.c_str()).getSpace().getSelectNpoints());
}

// ------------------------------------ read shape of the data -------------------------------------

inline std::vector<size_t> File::shape(std::string path)
{
  // check existence of path
  if ( ! exists(path) )
    throw std::runtime_error("HDF5pp::shape: dataset not found ('"+path+"')");

  // open data-space
  H5::DataSpace dataspace = m_file.openDataSet(path.c_str()).getSpace();

  // get the size in each direction
  // - read rank (a.k.a number of dimensions)
  int rank = dataspace.getSimpleExtentNdims();
  // - allocate as HDF5-type
  std::vector<hsize_t> dimsf(rank);
  // - allocate as vector
  std::vector<size_t> shape(rank);
  // - read
  dataspace.getSimpleExtentDims(dimsf.data(), NULL);
  // - convert to vector
  for ( int i = 0 ; i < rank ; ++i ) shape[i] = static_cast<size_t>(dimsf[i]);

  return shape;
}

// ------------------------- read shape of the data along a specific axis --------------------------

inline size_t File::shape(std::string path, size_t i)
{
  // check existence of path
  if ( ! exists(path) )
    throw std::runtime_error("HDF5pp::shape: dataset not found ('"+path+"')");

  // open dataset
  H5::DataSet   dataset   = m_file.openDataSet(path.c_str());
  H5::DataSpace dataspace = dataset.getSpace();

  // get the size in each direction
  // - read rank (a.k.a number of dimensions)
  int rank = dataspace.getSimpleExtentNdims();
  // - check rank
  if ( rank < static_cast<int>(i) )
    throw std::runtime_error("HDF5pp::shape: rank too low ('"+path+"')");
  // - allocate as HDF5-type
  std::vector<hsize_t> dimsf(rank);
  // - read
  dataspace.getSimpleExtentDims(dimsf.data(), NULL);

  return static_cast<size_t>(dimsf[i]);
}

// ----------------------- read the size of the data in an opened dataset ------------------------

inline size_t File::size(const H5::DataSet &dataset)
{
  return static_cast<size_t>(dataset.getSpace().getSelectNpoints());
}

// ----------------------- read the size of the data in an opened data-space -----------------------

inline size_t File::size(const H5::DataSpace &dataspace)
{
  return static_cast<size_t>(dataspace.getSelectNpoints());
}

// ------------------------ read the shape of the data in an opened dataset ------------------------

inline std::vector<size_t> File::shape(const H5::DataSet &dataset)
{
  // read the data-space
  H5::DataSpace dataspace = dataset.getSpace();

  // get the size in each direction
  // - read rank (a.k.a number of dimensions)
  int rank = dataspace.getSimpleExtentNdims();
  // - allocate as HDF5-type
  std::vector<hsize_t> dimsf(rank);
  // - allocate as vector
  std::vector<size_t> shape(rank);
  // - read
  dataspace.getSimpleExtentDims(dimsf.data(), NULL);
  // - convert to vector
  for ( int i = 0 ; i < rank ; ++i ) shape[i] = static_cast<size_t>(dimsf[i]);

  return shape;
}

// ---------------------- read the shape of the data in an opened data-space -----------------------

inline std::vector<size_t> File::shape(const H5::DataSpace &dataspace)
{
  // get the size in each direction
  // - read rank (a.k.a number of dimensions)
  int rank = dataspace.getSimpleExtentNdims();
  // - allocate as HDF5-type
  std::vector<hsize_t> dimsf(rank);
  // - allocate as vector
  std::vector<size_t> shape(rank);
  // - read
  dataspace.getSimpleExtentDims(dimsf.data(), NULL);
  // - convert to vector
  for ( int i = 0 ; i < rank ; ++i ) shape[i] = static_cast<size_t>(dimsf[i]);

  return shape;
}

// ============================= WRITE STD::STRING TO SEPARATE DATASET =============================

inline void File::write(std::string path, std::string input)
{
  // check existence of path
  if ( exists(path) )
    throw std::runtime_error("HDF5pp::write: path already exists ('"+path+"')");

  // create group(s) if needed
  createGroup(path);

  // set data-type and data-space
  H5::StrType   datatype(0, H5T_VARIABLE);
  H5::DataSpace dataspace(H5S_SCALAR);

  // create dataset
  H5::DataSet dataset = m_file.createDataSet(path.c_str(), datatype, dataspace);

  // write string to dataset
  dataset.write(input, datatype, dataspace);

  // flush the file if so requested
  if ( m_autoflush ) flush();
}

// =================== READ STD::STRING FROM DATASET THAT ONLY CONTAINS A STRING ===================

template<>
inline std::string File::read<std::string>(std::string path)
{
  // check existence of path
  if ( ! exists(path) )
    throw std::runtime_error("HDF5pp::read: dataset not found ('"+path+"')");

  // open dataset, get data-type
  H5::DataSet dataset = m_file.openDataSet(path.c_str());

  // allocate output
  std::string out;

  // read output
  dataset.read(out, dataset.getStrType(), dataset.getSpace());

  return out;
}

// ================= CHECK IF DATASET HAS A PRECISION THAT MATCHES A SPECIFIC TYPE =================

// ---------------------------------------------- int ----------------------------------------------

template<>
inline bool File::correct_presision<int>(const H5::DataSet &dataset)
{
  // check data-type
  if ( dataset.getTypeClass() != H5T_INTEGER ) return false;

  // check the number of bytes
  if ( dataset.getIntType().getSize() != sizeof(int) ) return false;

  return true;
}

// -------------------------------------------- size_t ---------------------------------------------

template<>
inline bool File::correct_presision<size_t>(const H5::DataSet &dataset)
{
  // check data-type
  if ( dataset.getTypeClass() != H5T_INTEGER ) return false;

  // check the number of bytes
  if ( dataset.getIntType().getSize() != sizeof(size_t) ) return false;

  return true;
}

// --------------------------------------------- float ---------------------------------------------

template<>
inline bool File::correct_presision<float>(const H5::DataSet &dataset)
{
  // check data-type
  if ( dataset.getTypeClass() != H5T_FLOAT ) return false;

  // check the number of bytes
  if ( dataset.getFloatType().getSize() != sizeof(float) ) return false;

  return true;
}

// -------------------------------------------- double ---------------------------------------------

template<>
inline bool File::correct_presision<double>(const H5::DataSet &dataset)
{
  // check data-type
  if ( dataset.getTypeClass() != H5T_FLOAT ) return false;

  // check the number of bytes
  if ( dataset.getFloatType().getSize() != sizeof(double) ) return false;

  return true;
}

// ================================ WRITE SCALAR TO SCALAR DATASET =================================

// ------------------------------------------- template --------------------------------------------

template<typename T>
inline void File::write(std::string path, T input, const H5::PredType& HT)
{
  // check existence of path
  if ( exists(path) )
    throw std::runtime_error("HDF5pp::write: path already exists ('"+path+"')");

  // create group(s) if needed
  createGroup(path);

  // define data-type, force little-endian storage
  auto datatype(HT);
  datatype.setOrder(H5T_ORDER_LE);

  // add dataset to file
  H5::DataSet dataset = m_file.createDataSet(path.c_str(), datatype, H5::DataSpace(H5S_SCALAR));

  // store data
  dataset.write(&input, HT);

  // flush the file if so requested
  if ( m_autoflush ) flush();
}

// ---------------------------------------------- int ----------------------------------------------

inline void File::write(std::string path, int input)
{
  return write<int>(path,input,H5::PredType::NATIVE_INT);
}

// -------------------------------------------- size_t ---------------------------------------------

inline void File::write(std::string path, size_t input)
{
  return write<size_t>(path,input,H5::PredType::NATIVE_HSIZE);
}

// --------------------------------------------- float ---------------------------------------------

inline void File::write(std::string path, float input)
{
  return write<float>(path,input,H5::PredType::NATIVE_FLOAT);
}

// -------------------------------------------- double ---------------------------------------------

inline void File::write(std::string path, double input)
{
  return write<double>(path,input,H5::PredType::NATIVE_DOUBLE);
}

// ============================== OVERWRITE SCALAR TO SCALAR DATASET ===============================

// ------------------------------------------- template --------------------------------------------

template<typename T>
inline void File::overwrite(std::string path, T input, const H5::PredType& HT)
{
  // new dataset: write using normal function
  if ( ! exists(path) ) return write<T>(path,input,HT);

  // check precision
  #ifndef HDF5PP_NDEBUG_PRECISION
    if ( ! this->correct_presision<T>(m_file.openDataSet(path.c_str())) )
      throw std::runtime_error("HDF5pp::overwrite: precision inconsistent ('"+path+"')");
  #endif

  // check size
  if ( size(path) != 1 )
    throw std::runtime_error("HDF5pp::overwrite: dataset has a rank different than 1 ('"+path+"')");

  // open dataset
  H5::DataSet *dataset = new H5::DataSet(m_file.openDataSet(path.c_str()));

  // store data
  dataset->write(&input, HT);

  // clean-up
  delete dataset;

  // flush the file if so requested
  if ( m_autoflush ) flush();
}

// ---------------------------------------------- int ----------------------------------------------

inline void File::overwrite(std::string path, int input)
{
  return overwrite<int>(path,input,H5::PredType::NATIVE_INT);
}

// -------------------------------------------- size_t ---------------------------------------------

inline void File::overwrite(std::string path, size_t input)
{
  return overwrite<size_t>(path,input,H5::PredType::NATIVE_HSIZE);
}

// --------------------------------------------- float ---------------------------------------------

inline void File::overwrite(std::string path, float input)
{
  return overwrite<float>(path,input,H5::PredType::NATIVE_FLOAT);
}

// -------------------------------------------- double ---------------------------------------------

inline void File::overwrite(std::string path, double input)
{
  return overwrite<double>(path,input,H5::PredType::NATIVE_DOUBLE);
}

// ======================== READ SCALAR FROM DATASET (SCALAR, OR OF SIZE 1) ========================

// ------------------------------------------- template --------------------------------------------

template<typename T>
inline T File::read_scalar(std::string path, const H5::PredType& HT)
{
  // check existence of path
  if ( ! exists(path) )
    throw std::runtime_error("HDF5pp::read: dataset not found ('"+path+"')");

  // open dataset
  H5::DataSet dataset = m_file.openDataSet(path.c_str());

  // check precision
  #ifndef HDF5PP_NDEBUG_PRECISION
    if ( ! this->correct_presision<T>(dataset) )
      throw std::runtime_error("HDF5pp::read: precision inconsistent ('"+path+"')");
  #endif

  // check size
  if ( this->size(dataset) > 1 )
    throw std::runtime_error("HDF5pp::read: dataset has a rank different than 1 ('"+path+"')");

  // allocate output
  T out;

  // read output
  dataset.read(&out, HT);

  return out;
}

// ---------------------------------------------- int ----------------------------------------------

template<>
inline int File::read<int>(std::string path)
{
  return read_scalar<int>(path,H5::PredType::NATIVE_INT);
}

// -------------------------------------------- size_t ---------------------------------------------

template<>
inline size_t File::read<size_t>(std::string path)
{
  return read_scalar<size_t>(path,H5::PredType::NATIVE_HSIZE);
}

// --------------------------------------------- float ---------------------------------------------

template<>
inline float File::read<float>(std::string path)
{
  return read_scalar<float>(path,H5::PredType::NATIVE_FLOAT);
}

// -------------------------------------------- double ---------------------------------------------

template<>
inline double File::read<double>(std::string path)
{
  return read_scalar<double>(path,H5::PredType::NATIVE_DOUBLE);
}

// ========================= WRITE SCALAR TO EXTENDABLE DATASET OF RANK 1 ==========================

// ------------------------------------------- template --------------------------------------------

template<typename T>
inline void File::write(std::string path, T input, const H5::PredType& HT,
  size_t index, T fill_val, size_t chunk_size
)
{
  // initialize if needed
  // --------------------

  if ( !exists(path) )
  {
    // create group(s) if needed
    createGroup(path);

    // set rank
    int rank = 1;

    // set initial and maximum shape of the array, and initial offset
    // NB initially an array of size "index+1" will be written, filled with "fill_val" and "input"
    std::vector<hsize_t> shape(rank, index+1), max_shape(rank, H5S_UNLIMITED), offset(rank, 0);

    // define the data-space
    H5::DataSpace dataspace(rank, shape.data(), max_shape.data());

    // choose chunk size (chosen by the user, who knows what to expect)
    std::vector<hsize_t> chunk_shape(rank, static_cast<hsize_t>(chunk_size));

    // enable chunking
    H5::DSetCreatPropList chunk_param;
    chunk_param.setChunk(rank, chunk_shape.data());
    chunk_param.setFillValue(HT, &fill_val);

    // create new dataset
    H5::DataSet dataset = m_file.createDataSet(path.c_str(), HT, dataspace, chunk_param);

    // make sure that the dataset is at least what we expect
    dataset.extend(shape.data());

    // select the hyperslap
    H5::DataSpace fspace = dataset.getSpace();
    fspace.selectHyperslab(H5S_SELECT_SET, shape.data(), offset.data());

    // initial data, filled with "fill_val" and with "input" at "init_data[index]"
    std::vector<T> init_data(index+1, fill_val);
    init_data[index] = input;

    // write data to the hyperslap
    dataset.write(init_data.data(), HT, dataspace, fspace);

    // flush the file if so requested
    if ( m_autoflush ) flush();

    // quit function
    return;
  }

  // extend
  // ------

  // open dataset
  H5::DataSet   dataset   = m_file.openDataSet(path.c_str());
  H5::DataSpace dataspace = dataset.getSpace();

  // check precision
  #ifndef HDF5PP_NDEBUG_PRECISION
    if ( ! this->correct_presision<T>(dataset) )
      throw std::runtime_error("HDF5pp::write: precision inconsistent ('"+path+"')");
  #endif

  // get the current rank and shape
  // - read rank (a.k.a number of dimensions)
  int rank = dataspace.getSimpleExtentNdims();
  // - allocate shape as HDF5-type
  std::vector<hsize_t> shape(rank);
  // - read shape
  dataspace.getSimpleExtentDims(shape.data(), NULL);

  // check rank (here only simple arrays are supported)
  if ( rank != 1 )
    throw std::runtime_error("HDF5pp::write: can only extend rank 1 array ('"+path+"')");

  // extend shape, if needed
  shape[0] = std::max(static_cast<size_t>(shape[0]), index+1);

  // process extension
  dataset.extend(shape.data());

  // set offset
  std::vector<hsize_t> offset(rank, index);

  // pseudo-shape of the array
  std::vector<hsize_t> extra_shape(rank, 1);

  // define the extra data-space
  H5::DataSpace extra_datasape(rank, extra_shape.data());

  // select a hyperslap
  H5::DataSpace fspace = dataset.getSpace();
  fspace.selectHyperslab(H5S_SELECT_SET, extra_shape.data(), offset.data());

  // write data to the hyperslap
  dataset.write(&input, HT, extra_datasape, fspace);

  // flush the file if so requested
  if ( m_autoflush ) flush();
}

// ---------------------------------------------- int ----------------------------------------------

inline void File::write(
  std::string path, int input, size_t index, int fill_val, size_t chunk_size
)
{
  return write<int>(path,input,H5::PredType::NATIVE_INT,index,fill_val,chunk_size);
}

// -------------------------------------------- size_t ---------------------------------------------

inline void File::write(
  std::string path, size_t input, size_t index, size_t fill_val, size_t chunk_size
)
{
  return write<size_t>(path,input,H5::PredType::NATIVE_HSIZE,index,fill_val,chunk_size);
}

// --------------------------------------------- float ---------------------------------------------

inline void File::write(
  std::string path, float input, size_t index, float fill_val, size_t chunk_size
)
{
  return write<float>(path,input,H5::PredType::NATIVE_FLOAT,index,fill_val,chunk_size);
}

// -------------------------------------------- double ---------------------------------------------

inline void File::write(
  std::string path, double input, size_t index, double fill_val, size_t chunk_size
)
{
  return write<double>(path,input,H5::PredType::NATIVE_DOUBLE,index,fill_val,chunk_size);
}

// ============================== READ SCALAR FROM DATASET OF RANK 1 ===============================

// ------------------------------------------- template --------------------------------------------

template<typename T>
inline T File::read(std::string path, const H5::PredType& HT, size_t index)
{
  // check existence of path
  if ( ! exists(path) )
    throw std::runtime_error("HDF5pp::read: dataset not found ('"+path+"')");

  // open dataset
  H5::DataSet   dataset   = m_file.openDataSet(path.c_str());
  H5::DataSpace dataspace = dataset.getSpace();

  // check precision
  #ifndef HDF5PP_NDEBUG_PRECISION
    if ( ! this->correct_presision<T>(dataset) )
      throw std::runtime_error("HDF5pp::read: incorrect precision ('"+path+"')");
  #endif

  // get the current rank and shape
  // - read rank (a.k.a number of dimensions)
  int rank = dataspace.getSimpleExtentNdims();
  // - allocate shape as HDF5-type
  std::vector<hsize_t> shape(rank);
  // - read shape
  dataspace.getSimpleExtentDims(shape.data(), NULL);

  // check rank (here only simple arrays are supported)
  if ( rank != 1 )
    throw std::runtime_error("HDF5pp::read: dataset not rank 1 ('"+path+"')");

  // check the shape
  if ( index >= shape[0] )
    throw std::runtime_error("HDF5pp::read: index out-of-bounds ('"+path+"')");

  // set offset
  std::vector<hsize_t> offset(rank, index);

  // pseudo-shape of the array
  std::vector<hsize_t> extra_shape(rank, 1);

  // define the extra data-space
  H5::DataSpace extra_datasape(rank, extra_shape.data());

  // select a hyperslap
  H5::DataSpace fspace = dataset.getSpace();
  fspace.selectHyperslab(H5S_SELECT_SET, extra_shape.data(), offset.data());

  // allocate output
  T out;

  // read data from the hyperslap
  dataset.read(&out, HT, extra_datasape, fspace);

  return out;
}

// ---------------------------------------------- int ----------------------------------------------

template<>
inline int File::read<int>(std::string path, size_t index)
{
  return read<int>(path,H5::PredType::NATIVE_INT,index);
}

// -------------------------------------------- size_t ---------------------------------------------

template<>
inline size_t File::read<size_t>(std::string path, size_t index)
{
  return read<size_t>(path,H5::PredType::NATIVE_HSIZE,index);
}

// --------------------------------------------- float ---------------------------------------------

template<>
inline float File::read<float>(std::string path, size_t index)
{
  return read<float>(path,H5::PredType::NATIVE_FLOAT,index);
}

// -------------------------------------------- double ---------------------------------------------

template<>
inline double File::read<double>(std::string path, size_t index)
{
  return read<double>(path,H5::PredType::NATIVE_DOUBLE,index);
}

// ====================== TEMPLATE TO WRITE ARRAY OF ARBITRARY SHAPE OR RANK =======================

template<typename T>
inline void File::write(
  std::string path, const T *input, const H5::PredType& HT, const std::vector<size_t> &shape
)
{
  // check existence of path
  if ( exists(path) )
    throw std::runtime_error("HDF5pp::write: path already exists ('"+path+"')");

  // create group(s) if needed
  createGroup(path);

  // shape of the array
  // - get the rank of the array
  size_t rank = shape.size();
  // - allocate shape as HDF5-type
  std::vector<hsize_t> dimsf(rank);
  // - store shape in each direction
  for ( size_t i = 0 ; i < rank ; ++i )
    dimsf[i] = shape[i];

  // define data-type, force little-endian storage
  auto datatype(HT);
  datatype.setOrder(H5T_ORDER_LE);

  // define data-space
  H5::DataSpace dataspace(rank, dimsf.data());

  // add dataset to file
  H5::DataSet dataset = m_file.createDataSet(path.c_str(), datatype, dataspace);

  // store data
  dataset.write(input, HT);

  // flush the file if so requested
  if ( m_autoflush ) flush();
}

// ==================== TEMPLATE TO OVERWRITE ARRAY OF ARBITRARY SHAPE OR RANK =====================

template<typename T>
inline void File::overwrite(
  std::string path, const T *input, const H5::PredType& HT, const std::vector<size_t> &shape
)
{
  // new dataset: write using normal function
  if ( ! exists(path) ) return write<T>(path,input,HT,shape);

  // check precision
  #ifndef HDF5PP_NDEBUG_PRECISION
    if ( ! this->correct_presision<T>(m_file.openDataSet(path.c_str())) )
      throw std::runtime_error("HDF5pp::overwrite: precision inconsistent ('"+path+"')");
  #endif

  // check shape
  if ( this->shape(path) != shape )
    throw std::runtime_error("HDF5pp::overwrite: shape inconsistent ('"+path+"')");

  // open dataset
  H5::DataSet *dataset = new H5::DataSet(m_file.openDataSet(path.c_str()));

  // store data
  dataset->write(input, HT);

  // clean-up
  delete dataset;

  // flush the file if so requested
  if ( m_autoflush ) flush();
}

// ======================= WRITE STD::VECTOR TO DATASET (OF ARBITRARY RANK) ========================

// ------------------------------------------- template --------------------------------------------

template<typename T>
inline void File::write(std::string path, const std::vector<T> &input, const H5::PredType& HT,
  const std::vector<size_t> &shape)
{
  // copy input shape
  std::vector<size_t> dims = shape;

  // default shape == size of input
  if ( dims.size() == 0 )
  {
    dims.resize(1);
    dims[0] = input.size();
  }

  // write to file
  write(path,input.data(),HT,dims);
}

// ---------------------------------------------- int ----------------------------------------------

inline void File::write(
  std::string path, const std::vector<int> &input, const std::vector<size_t> &shape
)
{
  return write(path,input,H5::PredType::NATIVE_INT,shape);
}

// -------------------------------------------- size_t ---------------------------------------------

inline void File::write(
  std::string path, const std::vector<size_t> &input, const std::vector<size_t> &shape
)
{
  return write(path,input,H5::PredType::NATIVE_HSIZE,shape);
}

// --------------------------------------------- float ---------------------------------------------

inline void File::write(
  std::string path, const std::vector<float> &input, const std::vector<size_t> &shape
)
{
  return write(path,input,H5::PredType::NATIVE_FLOAT,shape);
}

// -------------------------------------------- double ---------------------------------------------

inline void File::write(
  std::string path, const std::vector<double> &input, const std::vector<size_t> &shape
)
{
  return write(path,input,H5::PredType::NATIVE_DOUBLE,shape);
}

// ===================== OVERWRITE STD::VECTOR TO DATASET (OF ARBITRARY RANK) ======================

// ------------------------------------------- template --------------------------------------------

template<typename T>
inline void File::overwrite(std::string path, const std::vector<T> &input, const H5::PredType& HT,
  const std::vector<size_t> &shape)
{
  // copy input shape
  std::vector<size_t> dims = shape;

  // default shape == size of input
  if ( dims.size() == 0 )
  {
    dims.resize(1);
    dims[0] = input.size();
  }

  // write to file
  overwrite(path,input.data(),HT,dims);
}

// ---------------------------------------------- int ----------------------------------------------

inline void File::overwrite(
  std::string path, const std::vector<int> &input, const std::vector<size_t> &shape
)
{
  return overwrite(path,input,H5::PredType::NATIVE_INT,shape);
}

// -------------------------------------------- size_t ---------------------------------------------

inline void File::overwrite(
  std::string path, const std::vector<size_t> &input, const std::vector<size_t> &shape
)
{
  return overwrite(path,input,H5::PredType::NATIVE_HSIZE,shape);
}

// --------------------------------------------- float ---------------------------------------------

inline void File::overwrite(
  std::string path, const std::vector<float> &input, const std::vector<size_t> &shape
)
{
  return overwrite(path,input,H5::PredType::NATIVE_FLOAT,shape);
}

// -------------------------------------------- double ---------------------------------------------

inline void File::overwrite(
  std::string path, const std::vector<double> &input, const std::vector<size_t> &shape
)
{
  return overwrite(path,input,H5::PredType::NATIVE_DOUBLE,shape);
}

// ======================= READ STD::VECTOR FROM DATASET (OF ARBITRARY RANK) =======================

// ------------------------------------------- template --------------------------------------------

template<typename T>
inline std::vector<T> File::read_vector(std::string path, const H5::PredType& HT)
{
  // check existence of path
  if ( ! exists(path) )
    throw std::runtime_error("HDF5pp::read: dataset not found ('"+path+"')");

  // open dataset
  H5::DataSet dataset = m_file.openDataSet(path.c_str());

  // check precision
  #ifndef HDF5PP_NDEBUG_PRECISION
    if ( ! this->correct_presision<T>(dataset) )
      throw std::runtime_error("HDF5pp::read: precision inconsistent ('"+path+"')");
  #endif

  // allocate output
  std::vector<T> data(this->size(dataset));

  // read data
  dataset.read(const_cast<T*>(data.data()), HT);

  // return output
  return data;
}

// ---------------------------------------------- int ----------------------------------------------

template<>
inline std::vector<int> File::read<std::vector<int>>(std::string path)
{
  return read_vector<int>(path,H5::PredType::NATIVE_INT);
}

// -------------------------------------------- size_t ---------------------------------------------

template<>
inline std::vector<size_t> File::read<std::vector<size_t>>(std::string path)
{
  return read_vector<size_t>(path,H5::PredType::NATIVE_HSIZE);
}

// --------------------------------------------- float ---------------------------------------------

template<>
inline std::vector<float> File::read<std::vector<float>>(std::string path)
{
  return read_vector<float>(path,H5::PredType::NATIVE_FLOAT);
}

// -------------------------------------------- double ---------------------------------------------

template<>
inline std::vector<double> File::read<std::vector<double>>(std::string path)
{
  return read_vector<double>(path,H5::PredType::NATIVE_DOUBLE);
}

// ================================= WRITE EIGEN COLUMN TO DATASET =================================

#ifdef HDF5PP_EIGEN

// ------------------------------------------- template --------------------------------------------

template<typename T>
inline void File::write(std::string path,
  const Eigen::Matrix<T,Eigen::Dynamic,1,Eigen::ColMajor> &input, const H5::PredType& HT)
{
  // temporarily disable parallelization by Eigen (just in case)
  Eigen::setNbThreads(1);

  // set shape
  std::vector<size_t> shape(1, input.size());

  // write to file
  write(path,input.data(),HT,shape);

  // reset automatic parallelization by Eigen
  Eigen::setNbThreads(0);
}

// ---------------------------------------------- int ----------------------------------------------

inline void File::write(
  std::string path, const Eigen::Matrix<int,Eigen::Dynamic,1,Eigen::ColMajor> &input)
{
  return write(path,input,H5::PredType::NATIVE_INT);
}

// -------------------------------------------- size_t ---------------------------------------------

inline void File::write(
  std::string path, const Eigen::Matrix<size_t,Eigen::Dynamic,1,Eigen::ColMajor> &input)
{
  return write(path,input,H5::PredType::NATIVE_HSIZE);
}

// --------------------------------------------- float ---------------------------------------------

inline void File::write(
  std::string path, const Eigen::Matrix<float,Eigen::Dynamic,1,Eigen::ColMajor> &input)
{
  return write(path,input,H5::PredType::NATIVE_FLOAT);
}

// -------------------------------------------- double ---------------------------------------------

inline void File::write(
  std::string path, const Eigen::Matrix<double,Eigen::Dynamic,1,Eigen::ColMajor> &input)
{
  return write(path,input,H5::PredType::NATIVE_DOUBLE);
}

// -------------------------------------------------------------------------------------------------

#endif

// ================================= WRITE EIGEN MATRIX TO DATASET =================================

#ifdef HDF5PP_EIGEN

// ------------------------------------------- template --------------------------------------------

template<typename T>
inline void File::write(std::string path,
  const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &input,
  const H5::PredType& HT)
{
  // temporarily disable parallelization by Eigen (just in case)
  Eigen::setNbThreads(1);

  // set shape
  std::vector<size_t> shape(2);
  shape[0] = input.rows();
  shape[1] = input.cols();

  // write to file
  write(path,input.data(),HT,shape);

  // reset automatic parallelization by Eigen
  Eigen::setNbThreads(0);
}

// ---------------------------------------------- int ----------------------------------------------

inline void File::write(std::string path,
  const Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &input)
{
  return write(path,input,H5::PredType::NATIVE_INT);
}

// -------------------------------------------- size_t ---------------------------------------------

inline void File::write(std::string path,
  const Eigen::Matrix<size_t,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &input)
{
  return write(path,input,H5::PredType::NATIVE_HSIZE);
}

// --------------------------------------------- float ---------------------------------------------

inline void File::write(std::string path,
  const Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &input)
{
  return write(path,input,H5::PredType::NATIVE_FLOAT);
}

// -------------------------------------------- double ---------------------------------------------

inline void File::write(std::string path,
  const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &input)
{
  return write(path,input,H5::PredType::NATIVE_DOUBLE);
}

// -------------------------------------------------------------------------------------------------

#endif

// =============================== OVERWRITE EIGEN COLUMN TO DATASET ===============================

#ifdef HDF5PP_EIGEN

// ------------------------------------------- template --------------------------------------------

template<typename T>
inline void File::overwrite(std::string path,
  const Eigen::Matrix<T,Eigen::Dynamic,1,Eigen::ColMajor> &input, const H5::PredType& HT)
{
  // temporarily disable parallelization by Eigen (just in case)
  Eigen::setNbThreads(1);

  // set shape
  std::vector<size_t> shape(1, input.size());

  // overwrite to file
  overwrite(path,input.data(),HT,shape);

  // reset automatic parallelization by Eigen
  Eigen::setNbThreads(0);
}

// ---------------------------------------------- int ----------------------------------------------

inline void File::overwrite(
  std::string path, const Eigen::Matrix<int,Eigen::Dynamic,1,Eigen::ColMajor> &input)
{
  return overwrite(path,input,H5::PredType::NATIVE_INT);
}

// -------------------------------------------- size_t ---------------------------------------------

inline void File::overwrite(
  std::string path, const Eigen::Matrix<size_t,Eigen::Dynamic,1,Eigen::ColMajor> &input)
{
  return overwrite(path,input,H5::PredType::NATIVE_HSIZE);
}

// --------------------------------------------- float ---------------------------------------------

inline void File::overwrite(
  std::string path, const Eigen::Matrix<float,Eigen::Dynamic,1,Eigen::ColMajor> &input)
{
  return overwrite(path,input,H5::PredType::NATIVE_FLOAT);
}

// -------------------------------------------- double ---------------------------------------------

inline void File::overwrite(
  std::string path, const Eigen::Matrix<double,Eigen::Dynamic,1,Eigen::ColMajor> &input)
{
  return overwrite(path,input,H5::PredType::NATIVE_DOUBLE);
}

// -------------------------------------------------------------------------------------------------

#endif

// =============================== OVERWRITE EIGEN MATRIX TO DATASET ===============================

#ifdef HDF5PP_EIGEN

// ------------------------------------------- template --------------------------------------------

template<typename T>
inline void File::overwrite(std::string path,
  const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &input,
  const H5::PredType& HT)
{
  // temporarily disable parallelization by Eigen (just in case)
  Eigen::setNbThreads(1);

  // set shape
  std::vector<size_t> shape(2);
  shape[0] = input.rows();
  shape[1] = input.cols();

  // overwrite to file
  overwrite(path,input.data(),HT,shape);

  // reset automatic parallelization by Eigen
  Eigen::setNbThreads(0);
}

// ---------------------------------------------- int ----------------------------------------------

inline void File::overwrite(std::string path,
  const Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &input)
{
  return overwrite(path,input,H5::PredType::NATIVE_INT);
}

// -------------------------------------------- size_t ---------------------------------------------

inline void File::overwrite(std::string path,
  const Eigen::Matrix<size_t,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &input)
{
  return overwrite(path,input,H5::PredType::NATIVE_HSIZE);
}

// --------------------------------------------- float ---------------------------------------------

inline void File::overwrite(std::string path,
  const Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &input)
{
  return overwrite(path,input,H5::PredType::NATIVE_FLOAT);
}

// -------------------------------------------- double ---------------------------------------------

inline void File::overwrite(std::string path,
  const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &input)
{
  return overwrite(path,input,H5::PredType::NATIVE_DOUBLE);
}

// -------------------------------------------------------------------------------------------------

#endif

// ================================= READ TO DATASET EIGEN COLUMN =================================

#ifdef HDF5PP_EIGEN

// ------------------------------------------- template --------------------------------------------

template<typename T>
inline Eigen::Matrix<T,Eigen::Dynamic,1,Eigen::ColMajor> File::read_eigen_column(std::string path,
  const H5::PredType& HT)
{
  // check existence of path
  if ( ! exists(path) )
    throw std::runtime_error("HDF5pp::read: dataset not found ('"+path+"')");

  // open dataset
  H5::DataSet dataset = m_file.openDataSet(path.c_str());

  // check precision
  #ifndef HDF5PP_NDEBUG_PRECISION
    if ( ! this->correct_presision<T>(dataset) )
      throw std::runtime_error("HDF5pp::read: precision inconsistent ('"+path+"')");
  #endif

  // temporarily disable parallelization by Eigen (just in case)
  Eigen::setNbThreads(1);

  // allocate output
  Eigen::Matrix<T,Eigen::Dynamic,1,Eigen::ColMajor> data(this->size(dataset));

  // read data
  dataset.read(data.data(), HT);

  // reset automatic parallelization by Eigen
  Eigen::setNbThreads(0);

  // return output
  return data;
}

// ---------------------------------------------- int ----------------------------------------------

template<>
inline Eigen::Matrix<int,Eigen::Dynamic,1,Eigen::ColMajor>
File::read<Eigen::Matrix<int,Eigen::Dynamic,1,Eigen::ColMajor>>(std::string path)
{
  return read_eigen_column<int>(path,H5::PredType::NATIVE_INT);
}

// -------------------------------------------- size_t ---------------------------------------------

template<>
inline Eigen::Matrix<size_t,Eigen::Dynamic,1,Eigen::ColMajor>
File::read<Eigen::Matrix<size_t,Eigen::Dynamic,1,Eigen::ColMajor>>(std::string path)
{
  return read_eigen_column<size_t>(path,H5::PredType::NATIVE_HSIZE);
}

// --------------------------------------------- float ---------------------------------------------

template<>
inline Eigen::Matrix<float,Eigen::Dynamic,1,Eigen::ColMajor>
File::read<Eigen::Matrix<float,Eigen::Dynamic,1,Eigen::ColMajor>>(std::string path)
{
  return read_eigen_column<float>(path,H5::PredType::NATIVE_FLOAT);
}

// -------------------------------------------- double ---------------------------------------------

template<>
inline Eigen::Matrix<double,Eigen::Dynamic,1,Eigen::ColMajor>
File::read<Eigen::Matrix<double,Eigen::Dynamic,1,Eigen::ColMajor>>(std::string path)
{
  return read_eigen_column<double>(path,H5::PredType::NATIVE_DOUBLE);
}

// -------------------------------------------------------------------------------------------------

#endif

// ================================= READ TO DATASET EIGEN MATRIX ==================================

#ifdef HDF5PP_EIGEN

// ------------------------------------------- template --------------------------------------------

template<typename T>
inline Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>
File::read_eigen_matrix(std::string path, const H5::PredType& HT)
{
  // check existence of path
  if ( ! exists(path) )
    throw std::runtime_error("HDF5pp::read: dataset not found ('"+path+"')");

  // open dataset
  H5::DataSet dataset = m_file.openDataSet(path.c_str());

  // check precision
  #ifndef HDF5PP_NDEBUG_PRECISION
    if ( ! this->correct_presision<T>(dataset) )
      throw std::runtime_error("HDF5pp::read: precision inconsistent ('"+path+"')");
  #endif

  // get shape
  std::vector<size_t> shape = this->shape(dataset);

  // check rank
  if ( shape.size() != 2 )
    throw std::runtime_error("HDF5pp::read: dataset has a rank different than 2 ('"+path+"')");

  // temporarily disable parallelization by Eigen (just in case)
  Eigen::setNbThreads(1);

  // allocate output
  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> data(shape[0], shape[1]);

  // read data
  dataset.read(data.data(), HT);

  // reset automatic parallelization by Eigen
  Eigen::setNbThreads(0);

  // return output
  return data;
}

// ---------------------------------------------- int ----------------------------------------------

template<>
inline Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>
File::read<Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>(std::string path)
{
  return read_eigen_matrix<int>(path,H5::PredType::NATIVE_INT);
}

// -------------------------------------------- size_t ---------------------------------------------

template<>
inline Eigen::Matrix<size_t,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>
File::read<Eigen::Matrix<size_t,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>(std::string path)
{
  return read_eigen_matrix<size_t>(path,H5::PredType::NATIVE_HSIZE);
}

// --------------------------------------------- float ---------------------------------------------

template<>
inline Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>
File::read<Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>(std::string path)
{
  return read_eigen_matrix<float>(path,H5::PredType::NATIVE_FLOAT);
}

// -------------------------------------------- double ---------------------------------------------

template<>
inline Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>
File::read<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>(std::string path)
{
  return read_eigen_matrix<double>(path,H5::PredType::NATIVE_DOUBLE);
}

// -------------------------------------------------------------------------------------------------

#endif

// ======================= WRITE CPPMAT-ND-ARRAY TO DATASET OF MATCHING RANK =======================

#ifdef HDF5PP_CPPMAT

// ------------------------------------------- template --------------------------------------------

template<typename T>
inline void File::write(std::string path, const cppmat::array<T> &input, const H5::PredType& HT)
{
  write(path,input.data(),HT,input.shape());
}

// ---------------------------------------------- int ----------------------------------------------

inline void File::write(std::string path, const cppmat::array<int> &input)
{
  return write(path,input,H5::PredType::NATIVE_INT);
}

// -------------------------------------------- size_t ---------------------------------------------

inline void File::write(std::string path, const cppmat::array<size_t> &input)
{
  return write(path,input,H5::PredType::NATIVE_HSIZE);
}

// --------------------------------------------- float ---------------------------------------------

inline void File::write(std::string path, const cppmat::array<float> &input)
{
  return write(path,input,H5::PredType::NATIVE_FLOAT);
}

// -------------------------------------------- double ---------------------------------------------

inline void File::write(std::string path, const cppmat::array<double> &input)
{
  return write(path,input,H5::PredType::NATIVE_DOUBLE);
}

// -------------------------------------------------------------------------------------------------

#endif

// ===================== OVERWRITE CPPMAT-ND-ARRAY TO DATASET OF MATCHING RANK =====================

#ifdef HDF5PP_CPPMAT

// ------------------------------------------- template --------------------------------------------

template<typename T>
inline void File::overwrite(std::string path, const cppmat::array<T> &input, const H5::PredType& HT)
{
  overwrite(path,input.data(),HT,input.shape());
}

// ---------------------------------------------- int ----------------------------------------------

inline void File::overwrite(std::string path, const cppmat::array<int> &input)
{
  return overwrite(path,input,H5::PredType::NATIVE_INT);
}

// -------------------------------------------- size_t ---------------------------------------------

inline void File::overwrite(std::string path, const cppmat::array<size_t> &input)
{
  return overwrite(path,input,H5::PredType::NATIVE_HSIZE);
}

// --------------------------------------------- float ---------------------------------------------

inline void File::overwrite(std::string path, const cppmat::array<float> &input)
{
  return overwrite(path,input,H5::PredType::NATIVE_FLOAT);
}

// -------------------------------------------- double ---------------------------------------------

inline void File::overwrite(std::string path, const cppmat::array<double> &input)
{
  return overwrite(path,input,H5::PredType::NATIVE_DOUBLE);
}

// -------------------------------------------------------------------------------------------------

#endif

// ================================= READ TO DATASET CPPMAT MATRIX =================================

#ifdef HDF5PP_CPPMAT

// ------------------------------------------- template --------------------------------------------

template<typename T>
inline cppmat::array<T> File::read_cppmat_array(std::string path, const H5::PredType& HT)
{
  // check existence of path
  if ( ! exists(path) )
    throw std::runtime_error("HDF5pp::read: dataset not found ('"+path+"')");

  // open dataset
  H5::DataSet dataset = m_file.openDataSet(path.c_str());

  // check precision
  #ifndef HDF5PP_NDEBUG_PRECISION
    if ( ! this->correct_presision<T>(dataset) )
      throw std::runtime_error("HDF5pp::read: precision inconsistent ('"+path+"')");
  #endif

  // allocate output
  cppmat::array<T> data(this->shape(dataset));

  // read data
  dataset.read(data.data(), HT);

  // return output
  return data;
}

// ---------------------------------------------- int ----------------------------------------------

template<>
inline cppmat::array<int> File::read<cppmat::array<int>>(std::string path)
{
  return read_cppmat_array<int>(path,H5::PredType::NATIVE_INT);
}

// -------------------------------------------- size_t ---------------------------------------------

template<>
inline cppmat::array<size_t> File::read<cppmat::array<size_t>>(std::string path)
{
  return read_cppmat_array<size_t>(path,H5::PredType::NATIVE_HSIZE);
}

// --------------------------------------------- float ---------------------------------------------

template<>
inline cppmat::array<float> File::read<cppmat::array<float>>(std::string path)
{
  return read_cppmat_array<float>(path,H5::PredType::NATIVE_FLOAT);
}

// -------------------------------------------- double ---------------------------------------------

template<>
inline cppmat::array<double> File::read<cppmat::array<double>>(std::string path)
{
  return read_cppmat_array<double>(path,H5::PredType::NATIVE_DOUBLE);
}

// -------------------------------------------------------------------------------------------------

#endif

// ========================================= WRITE XTENSOR =========================================

#ifdef HDF5PP_XTENSOR

template<class E>
inline void File::write(std::string path, const xt::xexpression<E> &data)
{
  auto&& d_data = xt::eval(data.derived_cast());

  std::vector<size_t> shape(d_data.shape().cbegin(), d_data.shape().cend());

  write(path, d_data.begin(), getType<typename E::value_type>(), shape);
}

#endif

// ======================================= OVERWRITE XTENSOR =======================================

#ifdef HDF5PP_XTENSOR

template<class E>
inline void File::overwrite(std::string path, const xt::xexpression<E> &data)
{
  auto&& d_data = xt::eval(data.derived_cast());

  std::vector<size_t> shape(d_data.shape().cbegin(), d_data.shape().cend());

  overwrite(path, d_data.begin(), getType<typename E::value_type>(), shape);
}

#endif

// ========================================= READ XTENSOR ==========================================

#ifdef HDF5PP_XTENSOR

// -------------------------------------------------------------------------------------------------

template<class T>
inline auto File::xread(const std::string& path)
{
  return xread_impl<xt::xarray<T>>(path, getType<T>());
}

// -------------------------------------------------------------------------------------------------

template<class T, size_t N>
inline auto File::xread(const std::string& path)
{
  return xread_impl<xt::xtensor<T,N>>(path, getType<T>());
}

// -------------------------------------------------------------------------------------------------

template <class T>
inline T File::xread_impl(const std::string& path, const H5::PredType& HT)
{
  // check existence of path
  if ( ! exists(path) )
    throw std::runtime_error("HDF5pp::read: dataset not found ('"+path+"')");

  // open dataset
  H5::DataSet dataset = m_file.openDataSet(path.c_str());

  // allocate output
  T data = T::from_shape(this->shape(dataset));

  // read data
  dataset.read(data.begin(), HT);

  // return output
  return data;
}

// -------------------------------------------------------------------------------------------------

#endif

// =================================================================================================

} // namespace H5p

// -------------------------------------------------------------------------------------------------

#endif
