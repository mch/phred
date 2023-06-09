/* 
   phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004-2005 Matt Hughes <mhughe@uvic.ca>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#ifndef NETCDF_DATA_WRITER_H
#define NETCDF_DATA_WRITER_H

#include "DataWriter.hh"
#include "../Exceptions.hh"

#include "../config.h"

#ifdef USE_NETCDF
#include <netcdf.h>
#endif

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
 * NetCDF requires the time dimension to be the first dimension, and
 * the data varies fastest (is contiguous) along the last
 * dimension. This works out actually, because internally, grid data
 * varies fastest along z.
 *
 * \bug Needs to check for empty variable names, and make up names if
 * need be. Or throw an exception or something, since they need to be
 * unique names that can are identified in the Data object too. 
 *
 * \bug Need to check for existing variable names in the file. 
 */
class NetCDFDataWriter : public DataWriter {
private:
protected:
  string filename_;

  int ncid_; /**< NetCDF file id */
  int omode_; /**< File access mode to use. Defaults to
                 NC_WRITE|NC_SHARE */
  bool fopen_;

  // I should probably group all of this data into an object...

  typedef struct ncdfvar {
    string var_name_;
    vector<int> dim_ids_;
    vector<int> dim_lens_;
    int var_id_;
    bool time_dim_;
  } ncdfvar_t; 
  
  map<string, ncdfvar_t> vars_;

  bool clobber_; /**< Indicates if the file should be overwritten or
                    not. Defaults to false. */

  /**
   * Handle a NetCDF error... throws an exception. 
   */
  void handle_error(int status);

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

  /**
   * Implements the function from DataWriter. This calls the function
   * below. 
   */
  unsigned int write_data(unsigned int time_step, 
                          Variable &var, void *ptr, 
                          unsigned int len);  


  /**
   * Does recursive writing of packed data. This function will only
   * be called on rank 0. 
   */
  unsigned int write_data(int var_id, size_t *start, size_t *count, 
                          MPI_Datatype t, void *ptr, unsigned int len);  

public:

  /**
   * Create a new NetCDFWriter using the given filename. The file is
   * opened for writing. 
   *
   * @param filename the filename to use
   * @param clobber Overwrite the file when opening it. Default is false. 
   */
  NetCDFDataWriter(const char *filename, 
                   bool clobber = false);

  /**
   * Create a new NetCDFWriter
   */
  NetCDFDataWriter();

  ~NetCDFDataWriter();

  /**
   * Set the filename to write to. The file is opened for writing. 
   *
   * @param filename 
   */
  inline void set_filename(string filename)
  {
    filename_ = filename;
    init();
  }

  /**
   * Initialize this object; open the file. If the file is already
   * open, this simply returns. 
   */
  void init();

  /**
   * Deinit; close the file. If the file is closed, this simply
   * returns. 
   */
  void deinit();

  /**
   * Add a result that this data writer will have to know how to
   * handle. The file must be open (init() has to have been called)
   * before using this function. 
   *
   * @param result describes the result
   */
  void add_variable(Result &result);

  /**
   * Print a string representation to an ostream.
   */
  virtual ostream& to_string(ostream &os) const;

};
#else
class NetCDFDataWriter : public DataWriter {
private:

public:
  NetCDFDataWriter()
  {
    throw NoNetCDFException(); //("NetCDF Support is not available.");
  }

  NetCDFDataWriter(const char *filename, 
                   bool clobber = false)
  {
    throw NoNetCDFException(); //("NetCDF Support is not available.");
  }

  ~NetCDFDataWriter()
  {}

  inline void init()
  {}

  inline void deinit()
  {}

  inline void add_variable(Result &result)
  {}

  inline void set_filename(string filename)
  {}

  /**
   * Print a string representation to an ostream.
   */
  virtual ostream& to_string(ostream &os) const
  {
    return os << "NetCDFDataWriter: NetCDF support is not included in this build.";
  }

protected:
  virtual unsigned int write_data(unsigned int time_step, 
                                  Variable &var, MPI_Datatype t, 
                                  void *ptr, unsigned int len)
  {
    return 0;
  }

};
#endif

#endif // NETCDF_DATA_WRITER_H
