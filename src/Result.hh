#ifndef RESULT_H
#define RESULT_H

#include <string>
#include "Types.hh"
#include "Grid.hh"
#include "Data.hh"

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
 * is implied. 
 */
class Result
{
private:
  Result(const Result &rhs);
  const Result &operator=(const Result &rhs);

protected:
  string var_name_; /**< Variable name */
  Data data_;

  unsigned int num_dims_; /**< Number of dimensions */
  unsigned int *dim_lens_; /**< Dimension lengths */

  unsigned int time_start_; /**< Time step to start returning results at */
  unsigned int time_stop_; /**< Time step to stop returning results at */
  unsigned int time_space_; /**< Number of time steps to skip between
                               results. */ 

  string dw_name_; /**< DataWriter name we intend our results for */

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
    : num_dims_(0), dim_lens_(0), time_start_(0), time_stop_(~0),
      time_space_(0) //, dw_(0)
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
   * Looks at the grid and produces output
   *
   * @param grid a reference to a Grid object
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  virtual Data &get_result(Grid &grid, unsigned int time_step) = 0;

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
   * Called to perform any initialization that may be required. Does
   * nothing by default. 
   */
  virtual void init()
  {}

  /** 
   * Returns the number of dimensions
   * @return unsigned int, num dims
   */
  inline unsigned int get_num_dims()
  {
    return num_dims_;
  }

  /**
   * Returns the lengths of the dimensions
   * @return a pointer to the lengths of the dimensions, with the 
   * length returned by get_num_dims()
   */
  inline unsigned int *get_dim_lengths()
  {
    return dim_lens_;
  }
};

#endif // RESULT_H
