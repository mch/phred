#ifndef HWALL_H
#define HWALL_H

#include "BoundaryCondition.hh"

/**
 * Magnetic wall boundary condition. The H field components tangential
 * to the face are set to zero, but the one normal to the face is set
 * such that it's derivate at the face is zero. 
 */
class Hwall : public BoundaryCond
{
private:
protected:
  /**
   * Implements a loop across a GridPlane. 
   *
   * @param T a subclass of GridPlane
   * @param r the region_t to apply the condition to, usually found
   * using find_face()
   * @param grid the grid to apply the boundary condition to
   */
  template<class T>
  void condition(region_t r, Grid &grid);

public:
  Hwall();
  ~Hwall();

  /**
   * Apply the magnetic wall boundary condition to a face of a grid. 
   *
   * @param face the face to apply the boundary condition to. 
   * @param grid the grid to apply the boundary condition to. 
   * @param the field components to affect. 
   */
  void apply(Face face,
             Grid &grid, FieldType type);

};

#endif // HWALL_H
