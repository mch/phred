/** \class FreqGrid
 * \brief A grid which allows for frequency depended material properties. 
 *
 * A derivation of grid which allows for frequency dependent material,
 * including conductive materials which result in the imaginary part
 * of permittivity being dependent on frequency.
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
  field_t ***ex_sum_;
  field_t ***ey_sum_;
  field_t ***ez_sum_;

  /** 
   * Allocate memory for the grid. Called by setup_grid(). Calls
   * Grid::alloc_grid, but adds allocation for the running sums. 
   */ 
  virtual void alloc_grid();

  /**
   * Free the memory allocated for the grid. Calls Grid::free_grid(),
   * but also frees the memory allocated to keep track of the sums. 
   */ 
  virtual void free_grid();

public:
  FreqGrid();
  virtual ~FreqGrid();

  /**
   * Compute the next time step of the fields, taking into account the
   * frequency dispersive materials.
   */
  virtual void update_fields();

  /**
   * Calculate the material constants from the given material
   * library. Extra coefficients need to be calculated. 
   *
   * @param matlib the material library to load the materials from 
   */
  virtual void load_materials(MaterialLib &matlib);

  
};

#endif // FREQ_GRID_H
