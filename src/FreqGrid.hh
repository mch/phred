/** \class FreqGrid
 * \brief A grid which allows for frequency depended material properties. 
 *
 * A derivation of grid which allows for frequency dependent material,
 * including conductive materials which result in the imaginary part
 * of permittivity being dependent on frequency.
 *
 * \bug We'll do the "hard on memory" approach first, treating every
 * cell as if it were a plasma material, but with constants which
 * reduce the equations to the normal Taflove perfect conductor
 * equations. This will probably be faster than....
 *
 * \bug compute the normal update equations over the entire region
 * first and only re-compute the plasma equations in the region where
 * there is actually plasma... or...
 * 
 * \bug Put an "if" inside the loop to see if the plasma or normal
 * update equations should be used. This may kill pipeline
 * performance, but it might also be ok, because we can probably
 * expect to go for long stretches with the same "if" outcome. 
 */

#ifndef FREQ_GRID_H
#define FREQ_GRID_H

#include "Grid.hh"

class FreqGrid : public Grid
{
protected:
  /*
   * Precalculated contants for plasma materials. 0 for non-plamsa
   * materials, indexed in the same was as Ca_ etc.
   */

  mat_coef_t *vcdt_; /**< e ^ (Electron collision frequency times deltat_) */
  mat_coef_t *omegapsq_; /**< Plamsa frequency in rad/sec squared
                            times delta over vc */

  /*
   * Since only part of the region will be plasma, only store extra
   * plamsa bits for those regions and use a seperate index to
   * traverse them, which is ok since the loops always traverse the
   * grid in the same way.
   */

  field_t *dx_; /**< old ex for next time step */
  field_t *sx_; /**< Auxillary term */
  field_t *sxm1_; /**< Previous sx_ */
  field_t *sxm2_; /**< Previous sxm1_ */

  field_t *dy_; /**< old ey for next time step */
  field_t *sy_; /**< Auxillary term */
  field_t *sym1_; /**< Previous sy_ */
  field_t *sym2_; /**< Previous sym1_ */

  field_t *dz_; /**< old ez for next time step */
  field_t *sz_; /**< Auxillary term */
  field_t *szm1_; /**< Previous sz_ */
  field_t *szm2_; /**< Previous szm1_ */

  /**
   * Free the memory allocated for the grid. Calls Grid::free_grid(),
   * but also frees the memory allocated to keep track of the sums. 
   */ 
  virtual void free_grid();

  /**
   * Compute the update equatations for the Ex field component. 
   */
  virtual void update_ex();

  /**
   * Compute the update equatations for the Ey field component. 
   */
  virtual void update_ey();

  /**
   * Compute the update equatations for the Ez field component. 
   */
  virtual void update_ez();

public:
  FreqGrid();
  virtual ~FreqGrid();

  /**
   * Copy constructor. Pretty much everything is copied, but the grid
   * is condidered uninitialized. No memory is allocated, grid data
   * pointers are set to zero. 
   */
  FreqGrid(const FreqGrid &rhs);

  /**
   * Assignment operator. Pretty much everything is copied, but the
   * grid is condidered uninitialized. No memory is allocated, grid
   * data pointers are set to zero.
   */
  const FreqGrid &operator=(const FreqGrid &rhs);

  /**
   * Calculate the material constants from the given material
   * library. Extra coefficients need to be calculated. 
   *
   * @param matlib the material library to load the materials from 
   */
  virtual void load_materials(MaterialLib &matlib);

  /**
   * Deallocate the memory used to store material coeffcients and so
   * on. This function can only be used in define mode. 
   */
  virtual void free_material();
 
  /** 
   * Allocate memory for the grid. Called by setup_grid(). Calls
   * Grid::alloc_grid, but adds allocation for the running sums. 
   */ 
  virtual void alloc_grid();

};

#endif // FREQ_GRID_H
