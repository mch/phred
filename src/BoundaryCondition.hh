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

// Trust me, there is one, but #include "Grid.hh" is a circular include
class Grid;

/**
 * An abstract base class for boundary conditions. Subclass this to
 * implement real ones.
 */
class BoundaryCond
{
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
  region_t find_face(Face face, Grid &grid);

  /**
   * The thickness of the boundary condition. Usually zero, but will
   * be nonzero for PML's.
   */
  unsigned int thickness_;

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
