#ifndef PML_COMMON_H
#define PML_COMMON_H

#include "Types.hh"

/**
 * Data common to a set of PML's. This holds coefficients and stuff
 * that are used by all PML's. This is intended to be a member of a
 * Grid class, which will call the alloc_coeffs and free_coeffs
 * methods. Clients other the a Grid or a Pml won't be able to do much
 * with this...
 */
class PmlCommon 
{
  friend class Pml;
  friend class Grid;

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

  /**
   * Helper function to initialize ratios
   */
  void init_ratios(Face face, Grid &grid, Pml *p);

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
   * Set the common PML parameters and calculate coeffs and
   * stuff. Called by Grid when leaving define mode, when all the
   * boundary conditions have been set. 
   *
   * @param grid The grid the PML is being applied to. 
   */
  void init_coeffs(Grid &grid);
};

#endif
