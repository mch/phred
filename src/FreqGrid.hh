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
  // Running sums of field data, for dispersive or conductive material
  // Most materials of interest are non-magnetic, so we won't bother
  // with the H components for now.
  field_t *ex_sum_;
  field_t *ey_sum_;
  field_t *ez_sum_;

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

  /**
   * Compute the update equatations for the Hx field component. 
   */
  virtual void update_hx();

  /**
   * Compute the update equatations for the Hy field component. 
   */
  virtual void update_hy();

  /**
   * Compute the update equatations for the Hz field component. 
   */
  virtual void update_hz();

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
