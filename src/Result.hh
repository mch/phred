#ifndef RESULT_H
#define RESULT_H

#include <string>
#include <vector>

#include "Types.hh"
#include "Grid.hh"
#include "Data.hh"
#include "LifeCycle.hh"

using namespace std;

/**
 * An abstract base class used to implement results. This class
 * recieves a reference to a grid and a reference to a DataWriter. It
 * takes data from the grid, performs calculations, and then calls
 * functions in the DataWriter to save the data. 
 *
 * The method get_result() MUST be implemented. Is intended to be called
 * after every the updates for each time step have been computed. 
 *
 * get_result(), returns a structure containing an MPI derived
 * datatype describing the data, a pointer to the data, and the number
 * of elements in the result. If there is no data to be written, then
 * the number of elements is zero. 
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
  Result(const Result &rhs);
  const Result &operator=(const Result &rhs);

protected:
  string var_name_; /**< Variable name */
  Data data_;

  vector<int> dim_lens_; /**< Dimension lengths */
  vector<string> dim_names_; /**< Dimension names */

  unsigned int time_start_; /**< Time step to start returning results at */
  unsigned int time_stop_; /**< Time step to stop returning results at */
  unsigned int time_space_; /**< Number of time steps to skip between
                               results. */ 

  string dw_name_; /**< DataWriter name we intend our results for */

  bool time_dim_; /**< True if this result has a time dimension. If
                     false, DataWriters except only one result from
                     this Result, usually at time_stop_. */

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
    : var_name_("Result"), time_start_(0), time_stop_(~0),
      time_space_(0), time_dim_(true)//, dw_(0)
  {}

  virtual ~Result()
  {}

  /**
   * Returns true if this result has a time dimension, that it,
   * returns data at more than one time step. 
   */
  inline bool has_time_dimension()
  {
    return time_dim_;
  }

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
   * Looks at the grid and produces output
   *
   * @param grid a reference to a Grid object
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  virtual Data &get_result(const Grid &grid, unsigned int time_step) = 0;

  /**
   * Set the name of the variable. 
   * @param name a string with the name
   */
  inline void set_name(const string &name)
  {
    var_name_ = name;
    data_.set_var_name(name);
  }

  /**
   * Return the name of the variable. 
   * @return a string with the name in it. 
   */
  inline const string &get_name()
  {
    return var_name_;
  }

  /**
   * Returns the lengths of the dimensions
   * @return a reference to a vector of lengths of the dimensions
   */
  inline const vector<int> &get_dim_lengths()
  {
    return dim_lens_;
  }

  /**
   * Returns the names of the dimensions
   * @return a reference to the names of the dimensions
   */
  inline const vector<string> &get_dim_names()
  {
    return dim_names_;
  }

  /**
   * Subclasses can implement this to allocate memory or init
   * constants or whatever. Called just before the simulation
   * starts. The default implementation does nothing. 
   */
  //virtual void init(const Grid &grid)
  //{}
  
  /**
   * Subclasses can implement this to deallocate memory or
   * whatever. Called just after the simulation ends. The default
   * implementation does nothing. 
   */
  //virtual void deinit(const Grid &grid)
  //{}
};

#endif // RESULT_H
