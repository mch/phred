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

#ifndef UPML_BC_H
#define UPML_BC_H

#include "BoundaryCondition.hh"

/**
 * Uniaxial PML implementation, per Gedney. 
 */ 
class UPml : public BoundaryCond
{
private:
protected:
  /**
   * Stored flux values from the last time step, used to compute a
   * second-order accurate frequency dependence.  
   */ 
  field_t *d_;
  field_t *h_;

  // This a "shadow library" of material constants. These are
  // modified versions of the grid material library which add
  // conductivity, a little on the inner surface to quite a bit on
  // the outer surface. The first index is the thickness, the second
  // is the material number. 
  mat_coef_t **Ca_;
  mat_coef_t **Cbx_;
  mat_coef_t **Cby_;
  mat_coef_t **Cbz_;

  mat_coef_t **Da_;
  mat_coef_t **Dbx_;
  mat_coef_t **Dby_;
  mat_coef_t **Dbz_;


  mat_coef_t *sigmas_; /**< Changing sigma from the inside to
                          outside. */ 
  unsigned int poly_order_; /**< Order of the polynomial used to shape
                               the conductivity */


//   template<region_t region, field_t component, field_t curl1, field_t curl2>
//   void normal_update(Grid &grid);

//   template<region_t region, field_t component, field_t curl1, field_t curl2>
//   void upml_update(Grid &grid);

  void update_ex(Grid &grid);
  void update_ey(Grid &grid);
  void update_ez(Grid &grid);

  void update_hx(Grid &grid);
  void update_hy(Grid &grid);
  void update_hz(Grid &grid);

public:
  UPml();
  ~UPml();

  /**
   * Applies a boundary condition to a face of the grid.
   *
   * @param face the face to apply the boundary condition to. 
   * @param grid the grid to apply the boundary condition to. 
   * @param the field components to affect. 
   */
  virtual void apply(Face face,
                     Grid &grid, 
                     FieldType type);

  /**
   * Boundary condition type. Subclasses may implement. Only really
   * important for PMLs.
   */
  virtual BoundaryCondition get_type() 
  {
    return UPML;
  }

  /**
   * Init the boundary condition. 
   */ 
  virtual void init(const Grid &grid, Face face);

  /**
   * Cleanup the boundary condition. 
   */ 
  virtual void deinit(const Grid &grid, Face face);

  /**
   * Transfer 
   * 
   * @param sd the subdomain boundary condition that has to exchange
   * the data. 
   * @param bcface the face the PML is on
   * @param sdface the face the subdmoain is on
   */ 
  virtual void add_sd_bcs(SubdomainBc *sd, Face bcface, Face sdface);

};

#endif // UPML_BC_H
