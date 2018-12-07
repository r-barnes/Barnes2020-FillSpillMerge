#ifndef _array2d_hpp_
#define _array2d_hpp_

#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <netcdf.h>
#include <stdexcept>
#include <string>
#include <type_traits>

template<class T>
class Array2D {
 public:
  T *data;
 private:
  int mywidth  = -1;
  int myheight = -1;
  void getDimLength(const int ncid, const int dimnum);
  void loadDEM(const std::string filename);
  void loadNC(const std::string filename, const std::string datavar);
 public:
  Array2D() = default;
  Array2D(const std::string filename, const std::string datavar);
  Array2D(const int width0, const int height0, const T val0);
  ~Array2D();
  T&   operator()(const int x, const int y);
  T    operator()(const int x, const int y) const;
  T&   operator()(const int i);
  T    operator()(const int i) const;
  int  width()  const;
  int  height() const;
  int  size()   const;
  bool isEdgeCell(const int x, const int y) const;
  bool inGrid(const int x, const int y) const;
  void setAll(const T val0);
  int  xyToI(const int x, const int y) const;
};

template<class T>
void Array2D<T>::getDimLength(const int ncid, const int dimnum){
  char dimnamebuf[100];
  size_t dimlen;
  int retval;
  if((retval= nc_inq_dimname(ncid, dimnum, dimnamebuf)))
    throw std::runtime_error("Couldn't get name of dimension!");

  std::string dimname(dimnamebuf);

  if((retval=nc_inq_dimlen(ncid, dimnum, &dimlen)))
    throw std::runtime_error("Couldn't get length for dimension '"+dimname+"'!");

  if(dimname=="lat")
    myheight = dimlen;
  else if(dimname=="lon")
    mywidth = dimlen;
}

template<class T>
Array2D<T>::Array2D(const std::string filename, const std::string datavar){
  if(filename.substr(filename.size()-3)=="dem")
    loadDEM(filename);
  else if(filename.substr(filename.size()-2)=="nc")
    loadNC(filename, datavar);
  else
    throw std::runtime_error("Unrecognized file extension!");
}

template<class T>
void Array2D<T>::loadNC(const std::string filename, const std::string datavar){
  /* This will be the netCDF ID for the file and data variable. */
  int ncid, varid, retval, dim_count;

  /* Open the file. NC_NOWRITE tells netCDF we want read-only access
   * to the file.*/
  if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid)))
    throw std::runtime_error("Failed to open file '" + filename + "'!");

  if ((retval = nc_inq_ndims(ncid, &dim_count)))
    throw std::runtime_error("Failed to get number of dimensions from file '" + filename + "'!");

  if(dim_count!=3)
    throw std::runtime_error("File '" + filename + "' did not have 2 dimensions!");    

  getDimLength(ncid, 0);
  getDimLength(ncid, 1);
  getDimLength(ncid, 2);

  if(mywidth==-1 || myheight==-1)
    throw std::runtime_error("File '" + filename + "' did not have a lat or lon dimension!");    

  data = new float[mywidth*myheight];

  /* Get the varid of the data variable, based on its name. */
  if ((retval = nc_inq_varid(ncid, datavar.c_str(), &varid)))
    throw std::runtime_error("Failed to get dataset '"+datavar+"' from file '" + filename + "'!");

  /* Read the data. */
  //TODO: Check data type
  if ((retval = nc_get_var_float(ncid, varid, data)))
    throw std::runtime_error("Failed to read data from '"+datavar+" from file '" + filename + "'! Error: " + nc_strerror(retval));

 // for(int i=0;i<mywidth*myheight;i++)
   // std::cout<<data[i]<<std::endl;

  /* Close the file, freeing all resources. */
  if ((retval = nc_close(ncid)))
    throw std::runtime_error("Failed to close file '" + filename + "'!");
}



template<class T>
void Array2D<T>::loadDEM(const std::string filename){
  std::ifstream fin(filename);

  std::string header;
  int val;

  fin>>header>>mywidth;
  if(header!="ncols")
    throw std::runtime_error("Not ncols");

  fin>>header>>myheight;
  if(header!="nrows")
    throw std::runtime_error("Not ncols");

  fin>>header>>val; //xllcorner
  fin>>header>>val; //yllcorner
  fin>>header>>val; //cellsize
  fin>>header>>val; //no_data value

  data = new float[mywidth*myheight];

  for(int y=0;y<myheight;y++)
  for(int x=0;x<mywidth; x++)
    fin>>data[y*mywidth+x];
}




template<class T>
Array2D<T>::Array2D(const int width0, const int height0, const T val0) {
  data = new T[width0*height0];
  mywidth  = width0;
  myheight = height0;
  for(int i=0;i<mywidth*myheight;i++)
    data[i] = val0;
}

template<class T>
Array2D<T>::~Array2D() {
  delete[] data;
}


template<class T>
T& Array2D<T>::operator()(const int x, const int y) {
  return data[y*mywidth+x];
}

template<class T>
T  Array2D<T>::operator()(const int x, const int y) const {
  return data[y*mywidth+x];
}

template<class T>
T& Array2D<T>::operator()(const int i) {
  return data[i];
}

template<class T>
T  Array2D<T>::operator()(const int i) const {
  return data[i];
}

template<class T>
int Array2D<T>::width() const {
  return mywidth;
}

template<class T>
int Array2D<T>::height() const {
  return myheight;
}

template<class T>
int Array2D<T>::size() const {
  return mywidth*myheight;
}

template<class T>
bool Array2D<T>::isEdgeCell(const int x, const int y) const {
  return y==0 || y==myheight-1;
}

template<class T>
bool Array2D<T>::inGrid(const int x, const int y) const {
  return 0<=x && 0<=y && x<mywidth && y<myheight;
}

template<class T>
void Array2D<T>::setAll(const T val0){
  for(int i=0;i<size();i++)
    data[i] = val0;
}

template<class T>
int Array2D<T>::xyToI(const int x, const int y) const {
  assert(x>=0 && y>=0 && x<width() && y<height());
  return y*width()+x;
}













template<class T>
void SaveAsNetCDF(const Array2D<T> &arr, const std::string filename, const std::string datavar){
  /* When we create netCDF variables and dimensions, we get back an
  * ID for each one. */
  int ncid, x_dimid, y_dimid, varid;
  int dimids[2];

  //For error handling
  int retval;

  //Create the file. The NC_CLOBBER parameter tells netCDF to overwrite this
  //file, if it already exists.
  if ((retval = nc_create(filename.c_str(), NC_CLOBBER, &ncid)))
    throw std::runtime_error("Failed to create file '" + filename + "'! Error: " + nc_strerror(retval));

  //Define the dimensions. NetCDF will hand back an ID for each.
  if ((retval = nc_def_dim(ncid, "x", arr.height(), &x_dimid)))
    throw std::runtime_error("Failed to create x dimension in file '" + filename + "'! Error: " + nc_strerror(retval));
  if ((retval = nc_def_dim(ncid, "y", arr.width(), &y_dimid)))
    throw std::runtime_error("Failed to create y dimension in file '" + filename + "'! Error: " + nc_strerror(retval));

  //The dimids array is used to pass the IDs of the dimensions of the variable.
  dimids[0] = x_dimid;
  dimids[1] = y_dimid;

  //Define the variable. The type of the variable in this case is NC_INT (4-byte
  //integer).
  nc_type dtype;
  if(std::is_same<T, int32_t>::value)
    dtype = NC_INT;
  else if(std::is_same<T, float>::value)
    dtype = NC_FLOAT;
  else
    throw std::runtime_error("Unimplemented data type found when writing file '" + filename + "'!");

  if ((retval = nc_def_var(ncid, datavar.c_str(), dtype, 2, dimids, &varid)))
    throw std::runtime_error("Failed to create variable '" + datavar + "' in file '" + filename + "'! Error: " + nc_strerror(retval));

  //End define mode. This tells netCDF we are done defining metadata.
  if ((retval = nc_enddef(ncid)))
    throw std::runtime_error("Could not end define mode when making file '" + filename + "'! Error: " + nc_strerror(retval));

  //Write the pretend data to the file. Although netCDF supports reading and
  //writing subsets of data, in this case we write all the data in one
  //operation.
  if(std::is_same<T, int32_t>::value){
    if ((retval = nc_put_var_int(ncid, varid, (const int*)arr.data)))
      throw std::runtime_error("Failed to write data to file '" + filename + "'! Error: " + nc_strerror(retval));
  } else if(std::is_same<T, float>::value){
    if ((retval = nc_put_var_float(ncid, varid, (const float*)arr.data)))
      throw std::runtime_error("Failed to write data to file '" + filename + "'! Error: " + nc_strerror(retval));
  }

  //Close the file. This frees up any internal netCDF resources associated with
  //the file, and flushes any buffers.
  if ((retval = nc_close(ncid)))
    throw std::runtime_error("Failed to close file '" + filename + "'! Error: " + nc_strerror(retval));
}


#endif

