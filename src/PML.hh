#ifndef PML_H
#define PML_H

#include "BoundaryCondition.hh"

/**
 * Perfectly matched layers boundary conditions. Berenger did it
 * first, and then some other guys did it a different way. 
 *
 * \bug IMPLEMENT ME!
 */
class Pml : public BoundaryCond
{
private:
protected:
  
public:
  Pml();
  ~Pml();

  /**
   * Applys a PML boundary condition to a face of the grid. 
   *
   * @param face the face to apply the boundary condition to. 
   * @param grid the grid to apply the boundary condition to. 
   */
  void apply(Face face, Grid &grid);  
};

#endif // PML_H
