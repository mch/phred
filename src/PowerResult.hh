/* 
   phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004 Matt Hughes <mhughe@uvic.ca>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#ifndef POWER_RESULT_H
#define POWER_RESULT_H

#include "DFTResult.hh"
#include "GridPlane.hh"

/**
 * Calculates the power flowing through a surface at specific
 * frequencies. 
 *
 * \bug Only rectangular surfaces are currently supported.... results
 * should have tighter intergration with geometry objects so that it
 * would be possible to collect data on a specific face of a geometry
 * object.... geometries in general need to be re-thought. 
 */
class PowerResult : public DFTResult
{
public:

  PowerResult();
  PowerResult(field_t freq_start, field_t freq_stop, 
              unsigned int num_freqs);
  ~PowerResult();

  /**
   * Allocate memory, sanity checks. 
   */ 
  virtual void init(const Grid &grid);
  
  /**
   * Free memory etc
   */
  virtual void deinit(const Grid &grid);

  /**
   * Looks at the grid and produces output
   *
   * @param grid a reference to a Grid object
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  virtual map<string, Variable *> &get_result(const Grid &grid, 
                                              unsigned int time_step);

  /**
   * Set the plane to do the calculation on
   */
  inline void set_region(region_t r)
  {
    region_ = r;
  }

  

protected:
  field_t *freqs_; /**< Frequencies */ 
  field_t *power_real_; /**< Power at each frequency */ 
  field_t *power_imag_; /**< Power at each frequency */ 
  region_t region_; /**< The region to get the power through, should
                       be a plane, so the min and max on one axis
                       should be the same. */
  field_t cell_area_; /**< Area for each FDTD cell... */ 

  char step_x_; /**< The direction to step for H averaging (+1 or -1 only) */ 
  char step_y_; /**< The direction to step for H averaging (+1 or -1 only) */ 
  char step_z_; /**< The direction to step for H averaging (+1 or -1 only) */ 
  GridPlane *plane_; /**< The GridPlane we can use to access Grid data */
  Axis normal_; /**< The axis the plane is normal to */ 

  Variable real_var_; /**< Power data */
  Variable imag_var_; /**< Power data */
  Variable freq_var_; /**< Frequency data */

  field_t *et1r_;
  field_t *et1i_;
  field_t *ht1r_;
  field_t *ht1i_;

  field_t *et2r_;
  field_t *et2i_;
  field_t *ht2r_;
  field_t *ht2i_;

  unsigned int x_size_;
  unsigned int y_size_;
  unsigned int z_size_;

  inline unsigned int pi(unsigned int x, unsigned int y, 
                         unsigned int z) const
  {
    return z + (y + x*y_size_) * z_size_;
  }


private:
};

#endif
