/* 
   phred - Phred is a parallel finite difference time domain
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

#ifndef HDF5_DATA_WRITER_H
#define HDF5_DATA_WRITER_H

#include "DataWriter.hh"
#include "config.h"
#include "Exceptions.hh"

#ifdef USE_HDF5

#include <map>
#include <vector>

using namespace std;

// I don't want the class to disappear when using Python scripts, I
// just want it to throw an exception when you try to use it and you
// can't.  Then you can look for the exception in your script and try
// a different data writer.

/**
 * HDF5 data writer; collects all data to be written to one rank and
 * saves it to disk.
 *
 * \bug IMPLEMENT ME!
 */
class Hdf5DataWriter : public DataWriter {
private:
protected:
  string filename_;


public:
  Hdf5DataWriter(int rank, int size)
    : DataWriter(rank, size)
  {
    throw NoHDF5Exception(); 
  }

  ~Hdf5DataWriter()
  {}

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

class Hdf5DataWriter : public DataWriter {
public:
  Hdf5DataWriter(int rank, int size)
    : DataWriter(rank, size)
  {
    throw NoHDF5Exception(); 
  }

  ~Hdf5DataWriter()
  {}

  void init()
  {}

  void deinit()
  {}

  void add_variable(Result &result)
  {}

};

#endif // USE_HDF5

#endif // HDF5_DATA_WRITER_H
