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

#ifndef HWALL_H
#define HWALL_H

#include "BoundaryCondition.hh"

/**
 * Magnetic wall boundary condition. The H field components tangential
 * to the face are set to zero, but the one normal to the face is set
 * such that it's derivate at the face is zero. 
 */
class Hwall : public BoundaryCond
{
private:
protected:
//   /**
//    * Implements a loop across a GridPlane. 
//    *
//    * @param T a subclass of GridPlane
//    * @param r the region_t to apply the condition to, usually found
//    * using find_face()
//    * @param grid the grid to apply the boundary condition to
//    */
//   template<class T>
//   void condition(region_t r, Grid &grid);

public:
  Hwall();
  ~Hwall();

  /**
   * Apply the magnetic wall boundary condition to a face of a grid. 
   *
   * @param face the face to apply the boundary condition to. 
   * @param grid the grid to apply the boundary condition to. 
   * @param the field components to affect. 
   */
  void apply(Face face,
             Grid &grid, FieldType type);

  /**
   * Boundary condition type. Subclasses may implement. Only really
   * important for PMLs.
   */
  virtual BoundaryCondition get_type() const;

};

#endif // HWALL_H
