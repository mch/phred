#ifndef EWALL_H
#define EWALL_H

#include "BoundaryCondition.hh"

/**
 * Electric wall boundary condition. The E field components tangential
 * to the face are set to zero, but the one normal to the face is set
 * such that it's derivate at the face is zero. 
 */
class Ewall : public BoundaryCond
{
private:
protected:

  template<class T>
  void helper(region r, Grid &grid);

public:
  Ewall();
  virtual ~Ewall();

  /**
   * Apply the electric wall boundary condition to a face of a grid. 
   *
   * @param face the face to apply the boundary condition to. 
   * @param grid the grid to apply the boundary condition to. 
   */
  virtual void apply(Face face,
                     Grid &grid);

};

#endif // EWALL_H
