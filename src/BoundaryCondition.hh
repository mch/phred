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

#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include "Types.hh"
#include "config.h"
#include <assert.h>

// Trust me, there is one, but #include "Grid.hh" is a circular include
class Grid;
class SubdomainBc;

/**
 * An abstract base class for boundary conditions. Subclass this to
 * implement real ones.
 */
class BoundaryCond
{
  friend class Grid; // Gah. These shouldn't be friends! They should be ENEMIES!
private:

protected:
  /**
   * Returns the min and max coordinates along each axis for a given
   * face. Each boundary condition implements it's own loops, but
   * this factors out some of the tedious work in anycase. If the
   * thickness is nonzero, the returned region_t is thick.
   *
   * @param face the Face to work on, as defined in Types.hh
   * @param grid the Grid object that we are looking at
   * @return a region_t containing the coordinate mins and maxs
   */
  region_t find_face(Face face, const Grid &grid);

  /**
   * The thickness of the boundary condition. Usually zero, but will
   * be nonzero for PML's.
   */
  unsigned int thickness_;

  region_t bc_r_; /**< Boundary condition region in local coords */ 

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
    assert(x < bc_r_.xmax && y < bc_r_.ymax && z < bc_r_.zmax);
    //assert(z + (y + x*(bc_r_.ymax - bc_r_.ymin) 
    //            * (bc_r_.zmax - bc_r_.zmin)) < sz_);

    return z + (y + x*(bc_r_.ymax - bc_r_.ymin)) 
                * (bc_r_.zmax - bc_r_.zmin);
  }

  // These regions are regions *in the grid coordinates*, not in the
  // BC coordinates, which are the regions for each field component
  // to update. 

  region_t grid_r_; /**< Total Grid region */

  region_t grid_ex_r_; /**< Region over which ex update is applied. */
  region_t grid_ey_r_; /**< Region over which ey update is applied. */
  region_t grid_ez_r_; /**< Region over which ez update is applied. */

  region_t grid_hx_r_; /**< Region over which hx update is applied. */
  region_t grid_hy_r_; /**< Region over which hy update is applied. */
  region_t grid_hz_r_; /**< Region over which hz update is applied. */

  /**
   * Calculates the field component update regions. 
   */
  virtual void compute_regions(Face face, const Grid &grid);

  /** 
   * Compute the local boundary condition region to update for a
   * particular field component.
   */
  region_t find_local_region(region_t field_r);

public:
  BoundaryCond() : thickness_(0) {}
  virtual ~BoundaryCond() {}

  /**
   * Applies a boundary condition to a face of the grid.
   *
   * @param face the face to apply the boundary condition to. 
   * @param grid the grid to apply the boundary condition to. 
   * @param the field components to affect. 
   */
  virtual void apply(Face face,
                     Grid &grid, 
                     FieldType type) = 0;

  /**
   * Set the thickness of the PML. I.e. the number of cells devoted to
   * the PML along the normal to the face to which PML is being
   * applied, and into the grid. 
   *
   * @param thickness yup
   */
  void set_thickness(unsigned int thickness);

  /**
   * Returns the thickness of the boundary condition. 
   */
  unsigned int get_thickness();

  /**
   * Boundary condition type. Subclasses may implement. Only really
   * important for PMLs.
   */
  virtual BoundaryCondition get_type() 
  {
    return UNKNOWN;
  }

  /**
   * Our own special version of LifeCycle's init() which has an
   * additional parameter: the face number the bounary is on.
   */ 
  virtual void init(const Grid &grid, Face face)
  {}

  /**
   * Our own special version of LifeCycle's deinit() which has an
   * additional parameter: the face number the bounary is on.
   */ 
  virtual void deinit(const Grid &grid, Face face)
  {}

  /**
   * Boundaries may have to share data across subdomain
   * boundaries. This function is called when the grid is set up so
   * they can tell subdomain boundary conditions all about it.
   * 
   * @param sd the subdomain boundary condition that has to exchange
   * the data. 
   * @param bcface the face the PML is on
   * @param sdface the face the subdmoain is on
   */ 
  virtual void add_sd_bcs(SubdomainBc *sd, Face bcface, Face sdface)
  {}
};

/**
 * This class represents an unknown boundary condition. It's just here
 * so I can return references without worrying about dereferencing
 * null pointers in GridInfo.
 */
class UnknownBc : public BoundaryCond
{
public:
  UnknownBc()
  {}

  ~UnknownBc()
  {}

  void apply(Face face, Grid &grid, FieldType type)
  {}

};

#endif // BOUNDARY_CONDITION_H
