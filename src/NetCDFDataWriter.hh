#ifndef NETCDF_DATA_WRITER_H
#define NETCDF_DATA_WRITER_H

#include "DataWriter.hh"
#include <netcdf.h>

#include <exception>
#include <map>
#include <vector>

using namespace std;

/**
 * NetCDF data writer; collects all data to be written to one rank and
 * saves it to disk.
 */
class NetCDFDataWriter : public DataWriter {
private:
  NetCDFDataWriter();

protected:
  string filename_;

  int ncid_; /**< NetCDF file id */
  int omode_; /**< File access mode to use. Defaults to
                 NC_WRITE|NC_SHARE */
  bool fopen_;

  map<string, vector<int>> dim_ids_; /**< Maps variable names to the
                                        dimension id's in used by that
                                        variable. */

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
   */
  int get_dim(int i, int size);

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

  /**
   * Handle the data produced by a Result object. 
   *
   * @param data a Data object containing the data to handle
   */
  virtual void handle_data(unsigned int time_step, Data &data);

};

#endif // NETCDF_DATA_WRITER_H
