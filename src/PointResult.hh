#ifndef POINT_RESULT_H
#define POINT_RESULT_H

#include "Result.hh"
#include "Types.hh"

#include <mpi.h>

/** 
 * A simple class that outputs the value of the field components at a
 * specific point in space.
 */ 
class PointResult : public Result
{
private:
protected:
  /** 
   * The field data the we copy from the grid. In order: time, ex, ey, ez,
   * hx, hy, hz.
   */
  field_t field_data_[7];

  // Point in global space. Have to translate it to the local grid. 
  point_t point_;

public:

  PointResult();
  PointResult(point_t p);
  ~PointResult();

  /**
   * Look for the point we want to return, and return it. 
   *
   * @param grid a reference to a Grid object
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  virtual Data &get_result(Grid &grid, unsigned int time_step);

  /**
   * Set the point in global coordinates
   *
   * @param point
   */
  inline void set_point(point_t p)
  {
    point_ = p;
  }

  /**
   * Get the point in global coordinates
   *
   * @return point_t
   */
  inline point_t get_point()
  {
    return point_;
  }

};

#endif // POINT_RESULT_H
