#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include "Types.hh"

// Trust me, there is one, but #include "Grid.hh" is a circular include
class Grid;

/**
 * An abstract base class for boundary conditions. Subclass this to
 * implement real ones.
 */
class BoundaryCond
{
private:

protected:
  /**
   * Returns the min and max coordinates along each axis for a given
   * face. Each boundary condition implements it's own loops, but
   * this factors out some of the tedious work in anycase. 
   *
   * @param face the Face to work on, as defined in Types.hh
   * @param grid the Grid object that we are looking at
   * @return a region containing the coordinate mins and maxs
   */
  region find_face(Face face, Grid &grid);

  /**
   * The thickness of the boundary condition. Usually zero, but will
   * be nonzero for PML's.
   */
  unsigned int thickness_;

public:
  BoundaryCond() {}
  virtual ~BoundaryCond() {}

  /**
   * Applies a boundary condition to a face of the grid.
   *
   * @param face the face to apply the boundary condition to. 
   * @param grid the grid to apply the boundary condition to. 
   */
  virtual void apply(Face face,
                     Grid &grid) = 0;

  /**
   * Returns the thickness of the boundary condition. 
   */
  unsigned int get_thickness();

  /**
   * Creates a copy of the boundary condition object and returns a
   * pointer to it. The callers is responsible for freeing it later. 
   */
  //virtual BoundaryCond *clone() = 0;

};

#endif // BOUNDARY_CONDITION_H
