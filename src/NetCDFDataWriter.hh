#ifndef NETCDF_DATA_WRITER_H
#define NETCDF_DATA_WRITER_H

#include "DataWriter.hh"

#include "config.h"

#ifdef USE_NETCDF
#include <netcdf.h>
#endif

#include <exception>
#include <map>
#include <vector>

using namespace std;

// I don't want the class to disappear when using Python scripts, I
// just want it to throw an exception when you try to use it and you
// can't.  Then you can look for the exception in your script and try
// a different data writer.

#ifdef USE_NETCDF
/**
 * NetCDF data writer; collects all data to be written to one rank and
 * saves it to disk.
 *
 * \bug DANGEROUS! This should work for now, but it makes assumptions
 * that may not be all that valid for general MPI derived data types!
 */
class NetCDFDataWriter : public DataWriter {
private:
protected:
  string filename_;

  int ncid_; /**< NetCDF file id */
  int omode_; /**< File access mode to use. Defaults to
                 NC_WRITE|NC_SHARE */
  bool fopen_;

  map<string, vector<int> > dim_ids_; /**< Maps variable names to the
                                        dimension id's in used by that
                                        variable. */
  map<string, int> var_ids_; /**< Variable id's accessed by variable name **/ 

  /**
   * Handle a NetCDF error... throws an exception. 
   */
  virtual void handle_error(int status);

  /**
   * Get a dimension id. Dimensions are named by guessing that the
   * first dimension will be x, the second y, and the third z, and
   * additional dimensions will just be A. This function tests for
   * existing dimension names and lengths, and if it find ones that
   * matches, it is returned. 
   *
   * @param i dimension number
   * @param size dimension length
   * @param name dimension name, which is guessed from dimension
   * number if it is empty 
   */
  int get_dim(int i, int size, string name);

  virtual void *write_data(Data &data, MPI_Datatype t, 
                           void *ptr, unsigned int len);  


  /**
   * Does recursive writing of packed data. Increments the pointer
   * after writing each value, and returns it. This function will only
   * be called on rank 0. 
   */
  void *write_data(int var_id, size_t *start, size_t *count, 
                   MPI_Datatype t, void *ptr, unsigned int len);  

public:

  NetCDFDataWriter(int rank, int size);
  ~NetCDFDataWriter();

  /**
   * Set the filename to write to. 
   *
   * @param filename 
   */
  inline void set_filename(string filename)
  {
    filename_ = filename;
  }

  /**
   * Initialize this object; open the file. If the file is already
   * open, this simply returns. 
   */
  virtual void init();

  /**
   * Deinit; close the file. If the file is closed, this simply
   * returns. 
   */
  virtual void deinit();

  /**
   * Add a result that this data writer will have to know how to
   * handle. The file must be open (init() has to have been called)
   * before using this function. 
   *
   * @param result describes the result
   */
  virtual void add_variable(Result &result);
};
#else
class NetCDFDataWriter : public DataWriter {
public:
  NetCDFDataWriter(int rank, int size)
    : DataWriter(rank, size)
  {
    throw exception(); //("NetCDF Support is not available.");
  }

  ~NetCDFDataWriter()
  {}

  void init()
  {}

  void deinit()
  {}

  void add_variable(Result &result)
  {}
};
#endif

#endif // NETCDF_DATA_WRITER_H
