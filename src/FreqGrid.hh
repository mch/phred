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

/** \class FreqGrid
 * \brief A grid which allows for frequency depended material properties. 
 *
 * A derivation of grid which allows for frequency dependent material,
 * including conductive materials which result in the imaginary part
 * of permittivity being dependent on frequency.
 *
 * \bug We'll do the "hard on memory" approach first, treating every
 * cell as if it were a plasma material
 *
 * \bug Since it's unlikely for the entire region, or even a
 * significant porition of it to be plasma, a less memory-wastey
 * scheme should be used, but take care to update
 * setup_subdomain_data as well, since the MPI data transactions may
 * become difficult. 
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

  mat_coef_t *vcdt_; /**< e ^ (-Electron collision frequency times deltat_) */
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
  virtual void update_ex(region_t update_r);

  /**
   * Compute the update equatations for the Ey field component. 
   */
  virtual void update_ey(region_t update_r);

  /**
   * Compute the update equatations for the Ez field component. 
   */
  virtual void update_ez(region_t update_r);

  /**
   * This is called by set_define_mode(true), and it should add any
   * data that needs to be exchanged across subdomain boundaries to 
   * the SubdomainBc object. 
   *
   * @param sd the subdomain boundary to add data to be exchanged to.
   * @param face the face this subdomain boundary is on. 
   */
  virtual void setup_subdomain_data(SubdomainBc *sd, Face face);

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
  virtual void load_materials(shared_ptr<MaterialLib> matlib);

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
