#ifndef BLOCK_RESULT_H
#define BLOCK_RESULT_H

#include "Result.hh"
#include "Types.hh"

#include <mpi.h>

/**
 * This result produces a block of data within the grid. Generally it
 * will only return data for only one field component, but it can
 * optionally return them all. 
 */
class BlockResult : public Result
{
private:
protected:
  region_t region_;
  FieldComponent field_comp_;

  bool init_; /**< Set to true after init() has been called. */
public:
  BlockResult();
  BlockResult(region_t r, FieldComponent field_comp = FC_EY);
  ~BlockResult();

  /**
   * Return an MPI datatype with a pointer to the block of data. 
   *
   * @param grid a reference to a Grid object
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  Data &get_result(const Grid &grid, unsigned int time_step);

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

  /**
   * Set the field compoment to operate on
   *
   * @param field_comp field component
   */
  inline void set_field_component(FieldComponent field_comp)
  {
    field_comp_ = field_comp;
  }

  /**
   * Returns the field component we are returning
   *
   * @return a field component name
   */
  inline FieldComponent get_field_component()
  {
    return field_comp_;
  }

};

#endif // BLOCK_RESULT_H
