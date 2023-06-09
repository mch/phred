/* 
   Phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004-2005 Matt Hughes <mhughe@uvic.ca>

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

#ifndef PML_H
#define PML_H

#include "BoundaryCondition.hh"
#include "SubdomainBc.hh"
#include "PmlCommon.hh"
#include "../Data.hh"

#include <mpi.h>

#include "../config.h"
#include <assert.h>

/**
 * PML Variation types. Not that I know what they mean or anything. 
 */
enum PmlVariation_t {
  VC,
  VL,
  VP,
  VG
};

/**
 * Perfectly matched layers boundary conditions. Berenger did it
 * first, and then some other guys did it a different way. 
 *
 * if sigma/epsilon_0 = sigma_star/mu_0, then the wave impedance of
 * the lossy free-space medium equals that of lossless vacuum -> no
 * reflection. 
 *
 * This is based on Jan's implementation. I just reorganized the
 * code, factored out some common stuff, etc. 
 */
class Pml : public BoundaryCond
{
  friend class PmlCommon;
private:
protected:

  // PML Parameters (thickness is a member of BoundaryCond)

  PmlVariation_t variation_; /**< One of 'c', 'l', 'p', 'g' */
  float g_; /**< Variation type g has a float parameter */
  float nrml_refl_; /**< Normal reflection, must be 0 < nrml_refl_ < 100 */

  // PML coefficient intermediats; used by PmlCommon to compute the
  // actual coefficients
  float ratio_m_;
  float exponent_n_;
  float delta_bndy_;
  float geometric_delta_;
  int geometric_profile_;

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

  // MPI Derived data types for moving split field data across
  // subdomain boundaries
  MPI_Datatype z_vector_;
  MPI_Datatype y_vector_;
  MPI_Datatype x_vector_;
  MPI_Datatype xy_plane_;
  MPI_Datatype yz_plane_;
  MPI_Datatype xz_plane_;

  /** 
   * A helper function called by PmlCommon::init_ratios
   * @param x
   */
  float sigma_over_eps_int(float x);

  /**
   * Allocate memory for the field data. Thickness must be non zero. 
   */
  virtual void alloc_pml_fields(Face face, const Grid &grid);

  /**
   * Free the memory used to hold the PML field data.
   */
  virtual void free_pml_fields();

  /**
   * Implements the PML update equations on the face.
   */
  template<class T>
  void update(Grid &grid);

public:
  Pml();
  Pml(const Pml &rhs);
  Pml(PmlVariation_t variation, float g, float nrml_refl);
  Pml(PmlVariation_t variation, float nrml_refl);
  ~Pml();

  const Pml &operator=(const Pml &rhs);

  /**
   * Set up the PML; allocate storage space for fields, calculate
   * coefficients. 
   * @param face the face this PML is on
   * @param grid the grid this PML is on
   */
  virtual void init(const Grid &grid, Face face);

  /**
   * Applys a PML boundary condition to a face of the grid. 
   *
   * @param face the face to apply the boundary condition to. 
   * @param grid the grid to apply the boundary condition to. 
   * @param the field components to affect. 
   */
  virtual void apply(Face face, Grid &grid, FieldType type);  

  /**
   * Set the profile variation of the PML. May be one of 'c', 'l',
   * 'p', 'g'. If 'g', also call set_g_param() to set the required
   * parameter. 
   * @param variation 
   */
  inline void set_variation(PmlVariation_t variation)
  {
    variation_ = variation;
  }

  /**
   * Sets the parameter used by the 'g' profile variation type.
   * @param g the param
   */
  inline void set_g_param(float g)
  {
    g_ = g;
  }

  /** 
   * Set the normal reflection. Must be greater than 0.0 and less
   * than 100.0.
   * @param refl
   */
  inline void set_nrml_refl(float refl)
  {
    nrml_refl_ = refl;
  }

  /**
   * Adds RxTxData objects to the given subdomain boundary conditions
   * so that the split E and H field components can be exchanged when
   * necessary. 
   *
   * @param sd the subdomain boundary condition that has to exchange
   * the data. 
   * @param pmlface the face the PML is on
   * @param sdface the face the subdmoain is on
   */
  virtual void add_sd_bcs(SubdomainBc *sd, Face pmlface, Face sdface);

  /**
   * Returns a BoundayCondition type so the grid knows this is a pml.
   */
  virtual BoundaryCondition get_type() const;

protected: // Only called by apply().

  /**
   * Update the Ex field component inside the PML
   */
  void pml_update_ex(Grid &grid);

  /**
   * Update the Ey field component inside the PML
   */
  void pml_update_ey(Grid &grid);

  /**
   * Update the Ez field component inside the PML
   */
  void pml_update_ez(Grid &grid);

  /**
   * Update the Hx field component inside the PML
   */
  void pml_update_hx(Grid &grid);

  /**
   * Update the Hy field component inside the PML
   */
  void pml_update_hy(Grid &grid);

  /**
   * Update the Hz field component inside the PML
   */
  void pml_update_hz(Grid &grid);

};

#endif // PML_H
