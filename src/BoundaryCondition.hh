#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include "Grid.hh"

/**
 * An abstract base class for boundary conditions. Subclass this to
 * implement real ones.
 */
class BoundaryCondition
{
private:
protected:

public:
  BoundaryCondition() {}
  virtual ~BoundaryCondition() {}

  /**
   * Implement this method to apply the boundary condition. This is
   * kind of kludgy, but speed might be an issue. I didn't factor out
   * the loop like for the excitation, although that is still an
   * option. Doing that didn't seem to increase code size very much (2
   * out of 40 instructions), but that was a really trivial case, so
   * things may have gotten optimized out. 
   *
   * @param face the face to apply the boundary condition to. 
   * @param grid the grid to apply the boundary condition to. 
   */
  virtual void apply(Face face,
                     Grid &grid) = 0;

  /*
   * @param x_start starting point on x axis
   * @param x_end end point on x axis
   * @param y_start starting point on y axis
   * @param y_end end point on y axis
   * @param z_start starting point on z axis
   * @param z_end end point on z axis
   */
  /*
unsigned int x_start, unsigned int x_stop, 
                     unsigned int x_start, unsigned int x_stop, 
                     unsigned int x_start, unsigned int x_stop, 
  */
};

#endif // BOUNDARY_CONDITION_H
