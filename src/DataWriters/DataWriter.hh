/* 
   Phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004 Matt Hughes <mhughe@uvic.ca>

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

#ifndef DATA_WRITER_H
#define DATA_WRITER_H

#include <mpi.h>

#include <map>
#include <vector>

#include "../Data.hh"
#include "../Types.hh"
#include "../Results/Result.hh"
#include "../LifeCycle.hh"

using namespace std;

/**
 * An abstract base class that can be subclassed to create objects
 * which save data to disk, or perform other output
 * functions. DataWriter's have a list of variables, although some
 * writers may only be able to support one. Variables have dimensions,
 * units, and names. DataWriter's must support at least one
 * dimensional data.
 *
 * The subclass must implement add variable and provide it's own way
 * of handling per variable data. 
 */
class DataWriter : public LifeCycle
{
private:

protected:
  int rank_;  /**< MPI Rank of the process that handles data writing
                 (defaults to 0) */

  string filename_; /**< File to write data to */

  /** 
   * A helper function to gather data from all ranks onto rank 0 so
   * that it can write it out. DataWriters that use MPI-IO should not
   * use this function. If you do use this function by using the
   * default implementation of handle_data(), then override
   * write_data() and end_data() instead. 
   */ 
  void gather_data(unsigned int time_step, Variable &var);

   /**
    * This is implemented by subclasses to write a block of data in
    * whatever file format that particular DataWriter supports. This
    * function will only be called on rank 0. This does nothing by
    * default, but is defined because subclasses do not necessarially
    * have to override it (if they override handle_data() instead). 
    *
    * @param var The data object describing the data to write
    * @param t The MPI datatype to write; may be different than the one
    * in the Data object because we take the datatype apart
    * recursivly. It will be a simple MPI datatype, not a derived type. 
    * @param ptr A point to the data
    * @param len The number of bytes left to write
    *
    * @return the number of bytes written. 
    */
  virtual unsigned int write_data(unsigned int time_step, 
                                  Variable &var, MPI_Datatype t, 
                                  void *ptr, unsigned int len) = 0;

  /**
   * Returns a vector of MPI data types that describe how the data
   * coming from each node will fit into the buffer.
   *
   * @param var the variable object describing the stuff. 
   * @return a vector of mpi data types, one per node. 
   */ 
  vector<MPI_Datatype> gather_types(const Variable &var);

  /**
   * Returns the number of bytes that each node will be sending. 
   */ 
  vector<unsigned int> get_recieve_sizes(const Data &data);

public:
  DataWriter();

  virtual ~DataWriter();

  /**
   * Returns the rank that data is to be written out on. 
   */
  inline int get_rank()
  {
    return rank_;
  }


  /**
   * Set the filename to write to. The file is opened for writing. 
   *
   * @param filename 
   */
  virtual void set_filename(string filename);

  /**
   * Add a result that this data writer will have to know how to
   * handle. This may throw an exception if the DataWriter cannot
   * handle a particular Result for some reason (two results added
   * when only one is supported for example). 
   *
   * @param result describes the result
   */
  virtual void add_variable(Result &result) = 0;

  /**
   * Add a scalar to be written to the output file. This is intended
   * for small things that might be needed in post-processing.
   *
   * This throws a DataWriterException if the data writer does not
   * support it.
   *
   * @param name variable name in the output file.
   * @param value the value of the scalar
   */ 
  virtual void add_scalar(const char *name, double &value);

  /**
   * Handle the data produced by a Result object. The default
   * implementation marshalls data from all ranks to rank 0, and rank
   * 0 alone is responsible for writing the data to disk. 
   *
   * DataWriters which have the ability to use Parallel I/O should
   * override this function. 
   *
   * @param data a Data object containing the data to handle
   */
  virtual void handle_data(unsigned int time_step, 
                           map<string, Variable *> &data);

};

#endif // DATA_WRITER_H
