#ifndef DATA_WRITER_H
#define DATA_WRITER_H

#include "Types.hh"
#include <mpi.h>
#include "Data.hh"

#include <map>
#include <vector>

#include "Result.hh"
#include "LifeCycle.hh"

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
  DataWriter()
  {}

protected:
  int rank_;  /**< MPI Rank of this process */
  int size_; /**< Number of processes in MPI comm. */  

  /** 
   * A helper function to gather data from all ranks onto rank 0 so
   * that it can write it out. DataWriters that use MPI-IO should not
   * use this function. If you do use this function by using the
   * default implementation of handle_data(), then override
   * write_data() and end_data() instead. 
   */ 
  void gather_data(unsigned int time_step, Data &data);

   /**
   * Does recursive writing of packed data. This function will only
   * be called on rank 0. This does nothing by default, but is defined
   * because subclasses do not necessarially have to override it (if
   * they override handle_data() instead). 
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
  virtual unsigned int write_data(unsigned int time_step, 
                                  Data &data, MPI_Datatype t, 
                                  void *ptr, unsigned int len) = 0;
 
public:
  DataWriter(int rank, int size) 
    : rank_(rank), size_(size)
  { }

  virtual ~DataWriter()
  { }

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
   * Handle the data produced by a Result object. The default
   * implementation marshalls data from all ranks to rank 0, and rank
   * 0 alone is responsible for writing the data to disk. 
   *
   * DataWriters which have the ability to use Parallel I/O should
   * override this function. 
   *
   * @param data a Data object containing the data to handle
   */
  virtual void handle_data(unsigned int time_step, Data &data);

};

#endif // DATA_WRITER_H
