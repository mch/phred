/* 
   Phred - Phred is a parallel finite difference time domain
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

#ifndef RESULT_H
#define RESULT_H

#include <string>
#include <vector>

#include "../Exceptions.hh"
#include "../Types.hh"
#include "../Grid.hh"
#include "../Data.hh"
#include "../LifeCycle.hh"

using namespace std;

/**
 * Groups dimension information together. 
 */ 
typedef struct {
  unsigned int local_len_;
  unsigned int global_len_;
  string name_;
  unsigned int start_;
} Dimension;

/**
 * Contains all the data about a variable, including global and local
 * dimension lengths. Once a variable has been set up, information
 * describing the dimensions and type of the variable is considered
 * immutable. 
 *
 * A variable can choose to return no result however, by setting the
 * number of items at a given time step to zero. 
 *
 * \bug Now that we've got this pre and post variables stuff, we can
 * probably get rid of this has_time_dimension stuff, since stuff in
 * the variables_ map obviously is intended to have a time dimension.
 *
 * \bug This should probably be split into two classes inheriting from
 * a common base class. One type of variable should support common
 * multidimensional arrays, one should support structures (which this
 * doesn't really do that well right now, and which hasn't been used
 * at all), and perhaps there should be another for scalars. 
 */
class Variable 
{
  friend class Result;
  friend class DataWriter;
private:
protected:
  string var_name_; /**< Variable name */
  Data data_;

  vector<Dimension> dimensions_; /**< Dimension information for this
                                    variable. */ 

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
   * Add a dimension to this variable. Once set, dimensions are
   * considered immutable.
   *
   * @param name the name of the dimension, or some identifying string
   * @param local_length the *local* size of the dimention 
   * @param global_length the *global* size of the dimention 
   * @param start the starting point for the dimension in the global
   * scheme of things.   
   */
  inline void add_dimension(const char *name, unsigned int local_length, 
                            unsigned int global_length, 
                            unsigned int start)
  {
    Dimension new_dim;
    new_dim.local_len_ = local_length;
    new_dim.global_len_ = global_length;
    new_dim.name_ = name;
    new_dim.start_ = start;
    
    dimensions_.push_back(new_dim);
  }

  /**
   * Reset the variable. This should be called in init() by all
   * results, in case the variable was previously used.
   */ 
  inline void reset()
  {
    dimensions_.clear();
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
  inline MPI_Datatype get_element_type() const
  {
    return element_type_;
  }

  /**
   * Return a const reference to the data object
   */
  inline const Data &get_data() const
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
   * Returns the list of dimensions. 
   *
   * @return a const reference to a vector of dimensions. 
   */ 
  inline const vector<Dimension> &get_dimensions() const
  {
    return dimensions_;
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
  inline unsigned int get_num_ptrs() const 
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
  inline const void *get_ptr(unsigned int ptr_num = 0) const
  {
    return data_.get_ptr(ptr_num);
  }

  /**
   * Returns the number of items (MPI Datatype data blocks) contained
   * in each pointer. 
   */
  inline unsigned int get_num() const
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
  inline void set_ptr(const void *ptr, unsigned int ptr_num = 0)
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

  /**
   * A convenience function for quicker more concise variable set up. 
   *
   * @param name The name of the variable
   * param has_time True if this variable has a time dimension
   * @param datatype The MPI datatype used to access memory in the
   *                 sending process
   * @param ptr The pointer to the source data in the sending process
   * @param num The number of items of type datatype to send
   */ 
  inline void setup(const string &name, bool has_time, 
                    MPI_Datatype datatype, const void *ptr, 
                    unsigned int num)
  {
    set_name(name);
    has_time_dimension(has_time);
    set_datatype(datatype);
    set_ptr(ptr);
    set_num(num);
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
 */
class Result : public LifeCycle
{
private:

protected:
  // This map should hold shared_ptrs!
  map <string, Variable *> variables_; /**< Variables for this result. Most
                                          results will only have one
                                          variable. */  
  
  map <string, Variable *> pre_vars_; /**< Variables that are ready
                                         before time stepping
                                         start. */ 
  
  map <string, Variable *> post_vars_; /**< Variables that are ready
                                          after time stepping is
                                          done. */ 

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
        && (time_space_ == 0 
            || (time_step - time_start_) % time_space_ == 0))
      return true;
    else 
      return false; 
  }

  /**
   * Subclasses can implement this function to calculate results that
   * are made available before time stepping begins.
   */ 
  virtual void calculate_pre_result(const Grid &grid)
  {}

  /**
   * Subclasses can implement this function to calculate results that
   * need to be calculated at each time step. 
   */
  virtual void calculate_result(const Grid &grid, unsigned int time_step)
  {}

  /**
   * Subclasses can implement this function to calculate results that
   * can only be calculated once time stepping has been completed. 
   */ 
  virtual void calculate_post_result(const Grid &grid)
  {}


public:
  Result() 
    : time_start_(0), time_stop_(~0),
      time_space_(0)//, dw_(0)
  {}

  virtual ~Result()
  {}

  /**
   * Set the time related parameters
   *
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

  inline void set_time_start(unsigned int start)
  { time_start_ = start; }

  inline void set_time_stop(unsigned int stop)
  { time_stop_ = stop; }

  inline void set_time_space(unsigned int space)
  { time_space_ = space; }

  inline unsigned int get_time_start()
  { return time_start_; }

  inline unsigned int get_time_stop()
  { return time_stop_; }

  inline unsigned int get_time_space()
  { return time_space_; }


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
   * Returns the map of all variables, including pre, post, and time
   * variables. They are not expected to have any output available at
   * this point.
   */
  inline const map<string, Variable *> get_variables() const
  {
    map<string, Variable *> all;

    map<string, Variable *>::const_iterator iter, iter_e;

    iter = variables_.begin();
    iter_e = variables_.end();
    
    for (; iter != iter_e; ++iter)
      all[(*iter).first] = (*iter).second;

    iter = pre_vars_.begin();
    iter_e = pre_vars_.end();
    
    for (; iter != iter_e; ++iter)
      all[(*iter).first] = (*iter).second;

    iter = post_vars_.begin();
    iter_e = post_vars_.end();
    
    for (; iter != iter_e; ++iter)
      all[(*iter).first] = (*iter).second;

    return all;
  }

  /**
   * Looks at the grid and produces output. 
   *
   * @param grid a reference to a Grid object
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  map<string, Variable *> &get_result(const Grid &grid, 
                                      unsigned int time_step)
  { 
    calculate_result(grid, time_step);

    return variables_; 
  }
  
  /**
   * Pre-output. Returns a set of variables that are ready before time
   * stepping starts. Things like frequency and angles ranges, etc.
   *
   * @param grid the grid!
   */ 
  map<string, Variable *> &get_pre_result(const Grid &grid)
  { 
    calculate_pre_result(grid);

    return pre_vars_; 
  }

  /**
   * Post-output. Returns a set of variables that are ready after time
   * stepping ends. Things like far field transformations, etc. 
   *
   * @param grid the grid!
   */ 
  map<string, Variable *> &get_post_result(const Grid &grid)
  { 
    calculate_post_result(grid);

    return post_vars_; 
  }

  /**
   * Print a string representation to an ostream.
   */
  virtual ostream& to_string(ostream &os) const = 0;

  friend ostream& operator<< (ostream& os, const Result &r);
};

inline ostream& operator<< (ostream& os, const Result &r)
{
  return r.to_string(os);
}

#endif // RESULT_H
