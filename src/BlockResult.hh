#ifndef BLOCK_RESULT_H
#define BLOCK_RESULT_H

#include "Result.hh"
#include "Types.hh"

#include <mpi.h>

/**
 * This result produces a block of data within the grid. 
 */
class BlockResult : public Result
{
private:
protected:
  region_t region_;
  bool init_; /**< Set to true after init() has been called. */
public:
  BlockResult();
  BlockResult(region_t r);
  ~BlockResult();

  /**
   * Return an MPI datatype with a pointer to the block of data. 
   *
   * @param grid a reference to a Grid object
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  Data &get_result(Grid &grid, unsigned int time_step);

  /**
   * Set the region to return. 
   *
   * @param region
   */
  inline void set_region(region_t region)
  {
    region_ = region;
  }

  /**
   * Called to perform any initialization that may be required.
   * Converts the region to local grid coordinates and constructs the
   * MPI datatype.
   */
  virtual void init(const Grid &grid);

  /**
   * Returns the region this result deals with. 
   *
   * @return the region
   */
  inline region_t get_region()
  {
    return region_;
  }

};

#endif // BLOCK_RESULT_H
