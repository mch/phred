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

#ifndef RESULT_H
#define RESULT_H

#include <string>
#include <vector>

#include "Exceptions.hh"
#include "Types.hh"
#include "Grid.hh"
#include "Data.hh"
#include "LifeCycle.hh"

using namespace std;

/**
 * Contains all the data about a variable. 
 */
class Variable 
{
  friend class Result;
private:
protected:
  string var_name_; /**< Variable name */
  Data data_;

  vector<int> dim_lens_; /**< Dimension lengths */
  vector<string> dim_names_; /**< Dimension names */
  vector<unsigned int> dim_starts_; /**< Starting positions for the
                                       data in the final output. Used
                                       when the data is collected to
                                       one rank for writing, so it
                                       knows where the data from each
                                       rank goes. */

  bool time_dim_; /**< True if this variable has a time dimension. If
                     false, DataWriters except only one result from
                     this Result, usually at time_stop_. */

  MPI_Datatype element_type_; /**< The MPI data type of the individual
                                 elements of this variable. */ 
public:
  Variable()
    : var_name_("Result"), time_dim_(true), 
      element_type_(GRID_MPI_TYPE)
  {}

  ~Variable()
  {}

  /**
   * Add a dimension to this variable
   *
   * @param name the name of the dimension, or some identifying string
   * @param length the *local* size of the dimention 
   * @param start the starting point for the dimension in the global
   * scheme of things.   
   */
  inline void add_dimension(const char *name, unsigned int length, 
                            unsigned int start = 0)
  {
    dim_names_.push_back(name);
    dim_lens_.push_back(length);
    dim_starts_.push_back(start);
  }

  /**
   * Reset the variable. This should be called in init() by all
   * results, in case the variable was previously used.
   */ 
  inline void reset()
  {
    dim_names_.clear();
    dim_lens_.clear();
    dim_starts_.clear();
  }

  /**
   * Set the MPI data type for the individual elements. This MUST NOT
   * be a derived data type. 
   */ 
  inline void set_element_type(MPI_Datatype t)
  {
    if (t == MPI_CHAR || t == MPI_BYTE || t == MPI_SHORT || t == MPI_INT
        || t == MPI_LONG || t == MPI_FLOAT || t == MPI_DOUBLE
        || t == MPI_UNSIGNED_CHAR || t == MPI_UNSIGNED_SHORT
        || t == MPI_UNSIGNED || t == MPI_UNSIGNED_LONG
        || t == MPI_LONG_DOUBLE)
      element_type_ = t;
    else
      throw ResultException("Element data type must not be a MPI derived data type.");
  }

  /**
   * Returns the MPI data type for the individual elements of this
   * variable. 
   */
  inline MPI_Datatype get_element_type()
  {
    return element_type_;
  }

  /**
   * Return a const reference to the data object
   */
  inline const Data &get_data()
  {
    return data_;
  }

  /**
   * Returns true if this result has a time dimension, that it,
   * returns data at more than one time step. 
   */
  bool has_time_dimension() const
  {
    return time_dim_;
  }

  /**
   * Used to set wether or not a variable has a dimension that is time. 
   */
  void has_time_dimension(bool t)
  {
    time_dim_ = t;
  }  

  /**
   * Set the name of the variable. 
   * @param name a string with the name
   */
  inline void set_name(const string &name)
  {
    var_name_ = name;
  }

  /**
   * Return the name of the variable. 
   * @return a string with the name in it. 
   */
  inline const string &get_name() const
  {
    return var_name_;
  }

  /**
   * Returns the lengths of the dimensions
   * @return a reference to a vector of lengths of the dimensions
   */
  inline const vector<int> &get_dim_lengths() const
  {
    return dim_lens_;
  }

  /**
   * Returns the names of the dimensions
   * @return a reference to the names of the dimensions
   */
  inline const vector<string> &get_dim_names() const
  {
    return dim_names_;
  }

  /**
   * Returns the starting points of the dimensions for this rank. 
   */
  inline const vector<unsigned int> &get_dim_starts() const
  {
    return dim_starts_;
  }

  // The following functions are just helpers that forward to the
  // Data_ member...

  /**
   * Get the data type used to process memory. This will usually be a
   * derived data type.
   */
  inline MPI_Datatype get_datatype() const
  {
    return data_.get_datatype();
  }

  /**
   * Used to set the data type used to process memory. This will usually be a
   * derived data type.
   * @param type MPI derived data type
   */
  inline void set_datatype(MPI_Datatype type)
  {
    data_.set_datatype(type);
  }

  /**
   * Returns the number of pointers contained in this data
   * block. 
   */
  inline unsigned int get_num_ptrs()
  {
    return data_.get_num_ptrs();
  }

  /** 
   * Returns a pointer to some data. We make no guarantee about the
   * state of the pointer, because it is set by some other object. 
   *
   * @param ptr_num the pointer number to return, defaults to
   * zero. Must be less than the value returned by get_num_ptrs(). 
   */
  inline void *get_ptr(unsigned int ptr_num = 0)
  {
    return data_.get_ptr(ptr_num);
  }

  /**
   * Returns the number of items (MPI Datatype data blocks) contained
   * in each pointer. 
   */
  inline unsigned int get_num()
  {
    return data_.get_num();
  }

  /**
   * Sets a pointer. The optional parameter ptr_num determines which
   * pointer number to set. 
   *
   * @param ptr_num optional parameter the specifies the pointer
   * number to set. Defaults to zero. 
   */
  inline void set_ptr(void *ptr, unsigned int ptr_num = 0)
  {
    data_.set_ptr(ptr, ptr_num);
  }

  /**
   * Set the number of items. Setting to zero indicates no result. 
   */
  inline void set_num(unsigned int num)
  {
    data_.set_num(num);
  }
};

/**
 * An abstract base class used to implement results. This class
 * recieves a reference to a grid and a reference to a DataWriter. It
 * takes data from the grid, performs calculations, and then calls
 * functions in the DataWriter to save the data. 
 *
 * The method get_result() MUST be implemented. Is intended to be called
 * after every the updates for each time step have been computed. 
 *
 * get_result(), returns a vector of structures containing an MPI
 * derived datatype describing the data, a pointer to the data, and
 * the number of elements in the result. If there is no data to be
 * written, then the number of elements is zero. These structures
 * represent variables.
 *
 * The dimensions defined by the subclasses is NOT to include
 * time. Since this request gets called on once every time step, time
 * is implied. However, some results may accumulate data at each time
 * step and return only one chunck of data at the end, such as a
 * DFT. Such results can exclude themselves from the implicit time
 * dimension by setting time_dim = false in thier constructors. It's
 * up to DataWriters to hounor this flag of course. 
 *
 * If Results do exclude themselves from the time dimension, but
 * return data at more than one time step anyway, then later data
 * overwrites previous data, as far as the DataWriter is concerened. 
 *
 */
class Result : public LifeCycle
{
private:

protected:
  map <string, Variable *> variables_; /**< Variables for this result. Most
                                          results will only have one
                                          variable. */  

  unsigned int time_start_; /**< Time step to start returning results at */
  unsigned int time_stop_; /**< Time step to stop returning results at */
  unsigned int time_space_; /**< Number of time steps to skip between
                               results. */ 

  string dw_name_; /**< DataWriter name we intend our results for */
  
  string base_name_; /**< Base name to use for variables, usually set
                        by FDTD to match the name given to the result */ 

  /** 
   * Help subclasses know if they should return any results or not
   */
  inline bool result_time(unsigned int time_step) 
  {
    if (time_step >= time_start_ && time_step <= time_stop_
        && (time_space_ == 0 || (time_step - time_start_) % time_space_))
      return true;
    else 
      return false; 
  }

public:
  Result() 
    : time_start_(0), time_stop_(~0),
      time_space_(0)//, dw_(0)
  {}

  virtual ~Result()
  {}

  /**
   * Set the time related parameters
   * @param start time step to start returning results at
   * @param stop time step to stop returning results at
   * @param space number of time steps to skip between results
   */
  inline void set_time_param(unsigned int start, unsigned int stop, 
                             unsigned int space)
  {
    time_start_ = start;
    time_stop_ = stop;
    time_space_ = space; 
  }

  /**
   * Set the data writer name (intended where we read that from a
   * config file and set the actual data writer later). 
   */
  inline void set_dw_name(const string &dw)
  {
    dw_name_ = dw;
  }

  /**
   * Get the data writer name
   */
  inline const string &get_dw_name()
  {
    return dw_name_;
  }

  /**
   * Set the base name to use for variables generated by result. This
   * name will be concatenated to the start of the variable name to
   * make it unique in the output file.
   *
   * @param name
   */
  inline void set_name(const char *name)
  {
    base_name_ = name;
  }

  /**
   * Returns the string that is prefixed to variables generated by
   * this result.
   *
   * @return name
   */
  inline const string &get_name()
  {
    return base_name_;
  }

  /**
   * Returns the map of variables. They are not expected to have any
   * output available at this point.
   */
  inline const map<string, Variable *> &get_variables() const
  {
    return variables_;
  }

  /**
   * Looks at the grid and produces output
   *
   * @param grid a reference to a Grid object
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  virtual map<string, Variable *> &get_result(const Grid &grid, 
                                              unsigned int time_step) = 0;
};

#endif // RESULT_H
