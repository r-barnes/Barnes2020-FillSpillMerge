#ifndef _kr_netcdf_
#define _kr_netcdf_

#include <netcdf.h>
#include <richdem/common/Array2D.hpp>
#include <string>

namespace rd = richdem;

static void GetDimLength(const int ncid, const int dimnum, int &mywidth, int &myheight){
  char dimnamebuf[100];
  size_t dimlen;
  int retval;
  if((retval = nc_inq_dimname(ncid, dimnum, dimnamebuf)))
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
rd::Array2D<T> LoadNetCDF(const std::string filename, const std::string datavar){
  /* This will be the netCDF ID for the file and data variable. */
  int ncid, varid, retval, dim_count;

  /* Open the file. NC_NOWRITE tells netCDF we want read-only access
   * to the file.*/
  if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid))){
    std::cerr<<nc_strerror(retval)<<std::endl;
    throw std::runtime_error("Failed to open file '" + filename + "'!");
  }

  if ((retval = nc_inq_ndims(ncid, &dim_count)))
    throw std::runtime_error("Failed to get number of dimensions from file '" + filename + "'!");

  if(dim_count!=3)
    throw std::runtime_error("File '" + filename + "' did not have 2 dimensions!");    

  int mywidth;
  int myheight;
  GetDimLength(ncid, 0, mywidth, myheight);
  GetDimLength(ncid, 1, mywidth, myheight);
  GetDimLength(ncid, 2, mywidth, myheight);

  if(mywidth==-1 || myheight==-1)
    throw std::runtime_error("File '" + filename + "' did not have a lat or lon dimension!");    

  rd::Array2D<T> dem(mywidth,myheight);

  /* Get the varid of the data variable, based on its name. */
  if ((retval = nc_inq_varid(ncid, datavar.c_str(), &varid)))
    throw std::runtime_error("Failed to get dataset '"+datavar+"' from file '" + filename + "'!");

  /* Read the data. */
  //TODO: Check data type
  if ((retval = nc_get_var_float(ncid, varid, dem.data())))
    throw std::runtime_error("Failed to read data from '"+datavar+" from file '" + filename + "'! Error: " + nc_strerror(retval));

 // for(int i=0;i<mywidth*myheight;i++)
   // std::cout<<data[i]<<std::endl;

  /* Close the file, freeing all resources. */
  if ((retval = nc_close(ncid)))
    throw std::runtime_error("Failed to close file '" + filename + "'!");

  return dem;
}




template<class T>
rd::Array2D<T> LoadDEM(const std::string filename){
  std::ifstream fin(filename);

  if(!fin.good())
    throw std::runtime_error("Failed to open file '"+filename+"'!");

  std::string header;
  int val;

  int mywidth,myheight;

  fin>>header>>mywidth;
  if(header!="ncols")
    throw std::runtime_error("Not ncols");

  fin>>header>>myheight;
  if(header!="nrows")
    throw std::runtime_error("Not ncols");

  rd::Array2D<T> temp(mywidth,myheight);

  fin>>header>>val; //xllcorner
  fin>>header>>val; //yllcorner
  fin>>header>>val; //cellsize
  fin>>header>>val; //no_data value
  if(header=="NODATA_value")
    temp.setNoData(val);

  for(int y=0;y<myheight;y++)
  for(int x=0;x<mywidth; x++)
    fin>>temp(x,y);

  return temp;
}



template<class T>
rd::Array2D<T> LoadData(const std::string filename, const std::string datavar){
  if(filename.substr(filename.size()-3)=="dem")
    return LoadDEM<T>(filename);
  else if(filename.substr(filename.size()-2)=="nc")
    return LoadNetCDF<T>(filename, datavar);
  else
    throw std::runtime_error("Unrecognized file extension!");
}



template<class T>
void SaveAsNetCDF(const rd::Array2D<T> &arr, const std::string filename, const std::string datavar){
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
    if ((retval = nc_put_var_int(ncid, varid, (const int*)arr.data())))
      throw std::runtime_error("Failed to write data to file '" + filename + "'! Error: " + nc_strerror(retval));
  } else if(std::is_same<T, float>::value){
    if ((retval = nc_put_var_float(ncid, varid, (const float*)arr.data())))
      throw std::runtime_error("Failed to write data to file '" + filename + "'! Error: " + nc_strerror(retval));
  }

  //Close the file. This frees up any internal netCDF resources associated with
  //the file, and flushes any buffers.
  if ((retval = nc_close(ncid)))
    throw std::runtime_error("Failed to close file '" + filename + "'! Error: " + nc_strerror(retval));
}

#endif
