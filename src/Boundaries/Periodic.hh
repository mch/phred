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

#ifndef PERIODIC_H
#define PERIODIC_H

#include "BoundaryCondition.hh"
#include "../Excitations/PeriodicExcitation.hh"
#include "../Types.hh"

/**
 * This class implements periodic boundary conditions. An excitation
 * that will supply a plane wave must be given.
 *
 * That face to which this boundary condition is assigned always
 * recieves data, i.e. it copies data from the face on the opposite
 * side to itself.
 *
 * \bug Does not transfer data across MPI sub-domains!
 */ 
class Periodic : public BoundaryCond
{
public:

  /**
   * Construct a periodic bounadry object. This is a little different
   * from all the other BoundaryCond object. Create one object and set
   * that same object to all of the faces which will be periodic with
   * each other. The faces the boundary are applied to must be
   * perpendicular to the excitation plane. 
   *
   * @param signal the signal to excite in the excitation plane
   */ 
  Periodic(shared_ptr<PeriodicExcitation> pe);

  ~Periodic();

  /**
   * Applies a boundary condition to a face of the grid.
   *
   * @param face the face to apply the boundary condition to. 
   * @param grid the grid to apply the boundary condition to. 
   * @param the field components to affect. 
   */
  void apply(Face face, Grid &grid, FieldType type);
  
  /**
   * Boundary condition type. Subclasses may implement. Only really
   * important for PMLs.
   */
  inline BoundaryCondition get_type() const
  { return PERIODIC; }
  
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
  
private:

  shared_ptr<PeriodicExcitation> pe_; /**< Signal function to use as
                                         excitation */ 
  
  bool faces_[6]; /**< The faces to which this periodic boundary is
                     being applied. */ 

  bool valid_; /**< True if this periodic boundary set up is valid. */ 

  int exchange_rank_; /**< The MPI rank that we share data with (may
                         be the same as MPI_RANK) */ 

  MPI_Datatype exchange_type_; /**< The MPI datatype used to copy 
                                  data around. */

  void copy_e(Face face, Grid &grid);
  void copy_h(Face face, Grid &grid);
};

#endif // PERIODIC_H
