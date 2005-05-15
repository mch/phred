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

#ifndef ACSII_DATA_WRITER_H
#define ACSII_DATA_WRITER_H

#include "DataWriter.hh"
#include <fstream>
#include <exception>

using namespace std;

/**
 * Writes data to a text file that can be loaded by MATLAB. 
 */
class AsciiDataWriter : public DataWriter {
  
private:
protected:
  string filename_;
  string var_name_; /**< For future comparison. */

  ofstream file_;
  
  unsigned int dim_len_;
  bool time_dim_; /**< True if this variable has a time dimension */

  /**
   * Write an array of data
   * @param ptr location of start of array
   * @param len number of items in array
   * @return the advanced ptr
   */
  template<class T>
  unsigned int write_array(void *ptr, unsigned int len);

  /**
   * Does recursive writing of packed data. This function will only be
   * called on rank 0.
   */
  unsigned int write_data(MPI_Datatype t, void *ptr, 
                          unsigned int len);

  /**
   * Called by gather_data(), calls the above function. 
   */
  unsigned int write_data(unsigned int time_step, 
                          Variable &var, void *ptr, 
                          unsigned int len);

  /**
   * Called by gather_data() to indicate that all of the data
   * availble has been written. 
   */
  void end_data();

public:
  AsciiDataWriter();
  AsciiDataWriter(const char *fn, Result &result);
  ~AsciiDataWriter();

  /**
   * Set the filename
   *
   * @param filename
   */
  inline void set_filename(string filename)
  {
    filename_ = filename;
  }

  /**
   * Returns the filename
   *
   * @return a string with the filename
   */
  inline string get_filename(string filename)
  {
    return filename_;
  }

  /**
   * Initialize this object; open the file
   */
  void init(const Grid &grid);

  /**
   * Deinit; close the file.
   */
  void deinit();

  /**
   * Add a variable that we should know about. Throws an exception if
   * you try to add more than one. This call currently only supports
   * one dimensional data as well. 
   *
   * @param result describes the variable
   */
  void add_variable(Result &result);

  /**
   * Print a string representation to an ostream.
   */
  virtual ostream& to_string(ostream &os) const;

};

#endif // ACSII_DATA_WRITER_H
