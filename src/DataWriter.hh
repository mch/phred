#ifndef DATA_WRITER_H
#define DATA_WRITER_H

#include "Types.hh"
#include <mpi.h>
#include "Data.hh"

#include <map>
#include <vector>

#include "Result.hh"

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
class DataWriter 
{
private:
public:

protected:
  int rank_;  /**< MPI Rank of this process */
  int size_; /**< Number of processes in MPI comm. */

private:
  DataWriter()
  {}

public:
  DataWriter(int rank, int size) 
    : rank_(rank), size_(size) 
  {}

  virtual ~DataWriter()
  {}

  /**
   * Initialize this object
   */
  virtual void init() = 0;

  /**
   * Deinit
   */
  virtual void deinit() = 0;

  /**
   * Returns the rank
   */
  inline int get_rank()
  {
    return rank_;
  }

  /**
   * Returns the size
   *
   * @return an int; number of MPI processes
   */
  inline int get_size()
  {
    return size_;
  }

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
   * Handle the data produced by a Result object. 
   *
   * @param data a Data object containing the data to handle
   */
  virtual void handle_data(unsigned int time_step, Data &data) = 0;

};

#endif // DATA_WRITER_H
