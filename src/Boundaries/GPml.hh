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

#ifndef GPML_BOUNDARY_CONDITION_H
#define GPML_BOUNDARY_CONDITION_H

#include "BoundaryCondition.hh"

/**
 * Generalized PML, useful for bounding lossy media and absorbing
 * those pesky evanescent waves.
 */ 
class GPML : public BoundaryCond 
{
public:
  GPml();
  GPml(const GPml &rhs);
  ~GPml();

  const GPml &operator=(const GPml &rhs);

  /**
   * Applies a boundary condition to a face of the grid.
   *
   * @param face the face to apply the boundary condition to. 
   * @param grid the grid to apply the boundary condition to. 
   * @param the field components to affect. 
   */
  void apply(Face face, Grid &grid, FieldType type);

  /**
   * Our own special version of LifeCycle's init() which has an
   * additional parameter: the face number the bounary is on.
   */ 
  void init(const Grid &grid, Face face);

  /**
   * Our own special version of LifeCycle's deinit() which has an
   * additional parameter: the face number the bounary is on.
   */ 
  void deinit(const Grid &grid, Face face);

  /**
   * Boundary condition type.
   */
  inline BoundaryCondition get_type() const
  { return GPML; }

private:


};

#endif // GPML_BOUNDARY_CONDITION_H
