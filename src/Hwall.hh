#ifndef HWALL_H
#define HWALL_H

#include "BoundaryCondition.hh"

/**
 * Magnetic wall boundary condition. The H field components tangential
 * to the face are set to zero, but the one normal to the face is set
 * such that it's derivate at the face is zero. 
 */
class Hwall : public BoundaryCondition
{
private:
protected:

public:
  Hwall();
  virtual ~Hwall();

  /**
   * Apply the magnetic wall boundary condition to a face of a grid. 
   *
   * @param face the face to apply the boundary condition to. 
   * @param grid the grid to apply the boundary condition to. 
   */
  virtual void apply(Face face,
                     Grid &grid);

};

#endif // HWALL_H
