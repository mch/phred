#ifndef PML_H
#define PML_H

#include "BoundaryCondition.hh"
#include "PmlCommon.hh"

/**
 * PML Variation types. Not that I know what they mean or anything. 
 */
enum PmlVariation_t {
  C,
  L,
  P,
  G
};

/**
 * Perfectly matched layers boundary conditions. Berenger did it
 * first, and then some other guys did it a different way. 
 *
 * if sigma/epsilon_0 = sigma_star/mu_0, then the wave impedance of
 * the lossy free-space medium equals that of lossless vacuum -> no
 * reflection. 
 *
 * \bug IMPLEMENT ME!
 */
class Pml : public BoundaryCond
{
private:
protected:

  // Split field component data. The memory layout is the same as for
  // the grid.

  field_t *exy_;
  field_t *exz_;
  
  field_t *eyx_;
  field_t *eyz_;
  
  field_t *ezx_;
  field_t *ezy_;
  
  field_t *hxy_;
  field_t *hxz_;
  
  field_t *hyx_;
  field_t *hyz_;
  
  field_t *hzx_;
  field_t *hzy_;

  bool alloced_; // True when we've allocated our memory, false
                 // otherwise. Only intended to make calling
                 // alloc_fields multiple times consqeunce free.

  region_t pml_r_; // A region

  /**
   * Allocate memory for the field data. Thickness must be non zero. 
   */
  void alloc_pml_fields(Face face, Grid &grid);

  /**
   * Free the memory used to hold the PML field data.
   */
  void free_pml_fields();


  PmlVariation_t variation_;
  float normal_refl_;
  float PML_g_;

  /**
   * Implements the PML update equations on the face.
   */
  template<class T>
  void update(Grid &grid);

public:
  Pml();
  ~Pml();

  /**
   * Point Index: Calculate the index in the arrays of a 3d
   * coordinate. ALWAYS USE THIS FUNCTION, in case I change the way
   * things are organized for some reason. It's inline, so it should
   * compile out.
   *
   * @param x
   * @param y
   * @param z
   * @param an index into the field component and material arrays. 
   */
  inline unsigned int pi(unsigned int x, unsigned int y, 
                         unsigned int z)
  {
    return z + (y + x*(pml_r_.ymax - pml_r_.ymin) 
                * (pml_r_.zmax - pml_r_.zmin));
  }

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

  /**
   * Update the Ex field component inside the PML
   */
  void pml_update_ex(const region_t &grid_r, Grid &grid);

  /**
   * Update the Ey field component inside the PML
   */
  void pml_update_ey(const region_t &grid_r, Grid &grid);

  /**
   * Update the Ez field component inside the PML
   */
  void pml_update_ez(const region_t &grid_r, Grid &grid);

  /**
   * Update the Hx field component inside the PML
   */
  void pml_update_hx(const region_t &grid_r, Grid &grid);

  /**
   * Update the Hy field component inside the PML
   */
  void pml_update_hy(const region_t &grid_r, Grid &grid);

  /**
   * Update the Hz field component inside the PML
   */
  void pml_update_hz(const region_t &grid_r, Grid &grid);

};

#endif // PML_H
