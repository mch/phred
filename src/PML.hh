#ifndef PML_H
#define PML_H

#include "BoundaryCondition.hh"

// THIS IS A MESS RIGHT NOW. WORK STOPPED UNTIL I CAN GET A BETTER
// GRASP OF EXACTLY WHAT IS GOING ON WITH JAN's PML CODE

// 
/**
 * PML Variation types. Not that I know what they mean or anything. 
 */
enum PmlVariation_t {
  C,
  L,
  P,
  G
};

// forward decl
class Pml;

/**
 * Data common to a set of PML's. This holds coefficients and stuff
 * that are used by all PML's. 
 */
class PmlCommon {
  friend class Pml;

private:
  float *ratio_x_;
  float *ratio_star_x_;
  
  float *ratio_y_;
  float *ratio_star_y_;
  
  float *ratio_z_;
  float *ratio_star_z_;
  
  float *e_x_coef1_;
  float *e_x_coef2_;
  
  float *e_y_coef1_;
  float *e_y_coef2_;
  
  float *e_z_coef1_;
  float *e_z_coef2_;
  
  float *h_x_coef1_;
  float *h_x_coef2_;
  
  float *h_y_coef1_;
  float *h_y_coef2_;
  
  float *h_z_coef1_;
  float *h_z_coef2_;

  /**
   * Setup the coefficients 
   */
  void alloc_coeffs(Grid &grid);

  /**
   * Free the coefficients.
   */
  void free_coeffs();

public:
  /**
   * Constructor
   */
  PmlCommon();

  /**
   * Destructornator!
   */
  ~PmlCommon();
  
  /**
   * Set the common PML parameters and calculate coeffs and stuff
   */
  void init_coeffs(Face face, Pml &pml, Grid &grid);
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

  // Split field component data:
  field_t ***exy_;
  field_t ***exz_;
  
  field_t ***eyx_;
  field_t ***eyz_;
  
  field_t ***ezx_;
  field_t ***ezy_;
  
  field_t ***hxy_;
  field_t ***hxz_;
  
  field_t ***hyx_;
  field_t ***hyz_;
  
  field_t ***hzx_;
  field_t ***hzy_;
  
  /**
   * Allocate memory for the field data
   *
   * @param r a region_t describing the size of the PML. 
   */
  void alloc_pml_fields(region_t r);

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
