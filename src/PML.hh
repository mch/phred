#ifndef PML_H
#define PML_H

#include "BoundaryCondition.hh"

/**
 * Perfectly matched layers boundary conditions. Berenger did it
 * first, and then some other guys did it a different way. 
 *
 * if sigma/epsilon_0 = sigma_star/mu_0, then the wave impedance of
 * the lossy free-space medium equals that of lossless vacuum -> no
 * reflection. 
 *
 * This boundary condition implements three set's of update
 * equations; one for faces, where there is only one PML, one for
 * edges, where two PML's overlap, and one for corners, where three
 * PML's overlap. 
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

  /**
   * Set the thickness of the PML. I.e. the number of cells devoted to
   * the PML along the normal to the face to which PML is being
   * applied, and into the grid. 
   *
   * @param thickness yup
   */
  void set_thickness(unsigned int thickness);

};

#endif // PML_H
