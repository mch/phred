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

#ifndef PY_DATA_WRITER_H
#define PY_DATA_WRITER_H

#ifdef USE_PY_BINDINGS
#include <Python.h>
#endif

#include "DataWriter.hh"
#include "Exceptions.hh"

#include "config.h"

#include <map>
#include <vector>

using namespace std;

/**
 * This data writer makes data available to DataWriters written in
 * Python, as derived from this class.
 */
class PyDataWriter : public DataWriter {
private:
  inline void set_filename(string filename)
  {}

protected:
  vector<string> varnames_; /**< Variable names; data is placed into a
                               Python dictionary in this object
                               reference by this same name. */
  
public:
  PyDataWriter(int rank, int size);
  ~PyDataWriter();

#ifdef USE_PY_BINDINGS


  /**
   * Init the PY writer, open file and stuff. 
   */
  void init();

  /**
   * Clean up, but don't destruct. 
   */
  void deinit();

  /**
   * Add a result that this data writer will have to know how to
   * handle. Throws an exception if the variable is not a
   * BlockResult, since only those are handled for now. 
   *
   * @param result describes the result
   */
  void add_variable(Result &result);

  /**
   * Writes data to a PY XML file. Called by the default DataWriter
   * implementation of handle_data(). 
   *
   * @param data The data object describing the data to write
   * @param t The MPI datatype to write; may be different than the one
   * in the Data object because we take the datatype apart
   * recursivly. 
   * @param ptr A point to the data
   * @param len The number of bytes left to write
   *
   * @return the number of bytes written. 
   */
  unsigned int write_data(unsigned int time_step, 
                          Data &data, MPI_Datatype t, 
                          void *ptr, unsigned int len);
#else
  inline void init()
  {}

  inline void deinit()
  {}

  inline void add_variable(Result &result)
  {}

  inline unsigned int write_data(unsigned int time_step, 
                                 Data &data, MPI_Datatype t, 
                                 void *ptr, unsigned int len)
  {
    return 0;
  }
#endif  
};

#endif // PY_DATA_WRITER_H
