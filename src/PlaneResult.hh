#ifndef PLANE_RESULT_H
#define PLANE_RESULT_H

#include "Result.hh"
#include "Types.hh"

#include <mpi.h>

/** 
 * A simple class that outputs the value of the field components at a
 * specific plane in space.
 */ 
class PlaneResult : public Result
{
private:
protected:
  // Position of the plane in global space. Have to translate it to
  // the local grid.
  point_t plane_;

  // The face the plane is parallel to
  Face face_;

  // The field component we are interested in. 
  FieldComponent field_;

public:

  PlaneResult();
  ~PlaneResult();

  /**
   * Look for the plane we want to return, and return it. 
   *
   * @param grid a reference to a Grid object
   * @return a reference to a data object, which contains an MPI
   * derived data type, a planeer, and the number of items in the
   * result.
   */
  virtual Data &get_result(Grid &grid, unsigned int time_step);

  /**
   * Set the plane position in global coordinates
   *
   * @param p
   */
  inline void set_plane(point_t p, Face face)
  {
    plane_ = p;
    face_ = face;
  }

  /**
   * Get the position of the plane in global coordinates
   *
   * @return point_t
   */
  inline point_t get_plane()
  {
    return plane_;
  }

  /**
   * Returns the face the plane is perpendicular to. 
   * @return Face
   */
  inline Face get_face()
  {
    return face_;
  }

  /** 
   * Set the field component to return. PlaneResult returns Ey by default. 
   */
  inline void set_field(FieldComponent field)
  {
    field_ = field;
  }

  /**
   * Initalize the result. Set's the size of the plane from the grid.
   */
  virtual void init(const Grid &grid);
};

#endif // PLANE_RESULT_H
