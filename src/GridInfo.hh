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

#ifndef GRID_INFO_H
#define GRID_INFO_H

#include "Types.hh"
#include "Boundaries/BoundaryCondition.hh"
#include "Boundaries/Ewall.hh"
#include "Boundaries/Hwall.hh"
#include "Boundaries/SubdomainBc.hh"
#include "Boundaries/Pml.hh"

#include <sys/types.h>
#include <boost/shared_ptr.hpp>
using namespace boost;

/**
 * This class is plain old data that is in class form "just in
 * case". Think of it as a structure. All members are public, but
 * this may change in the future. 
 *
 * This class holds data about grids. Number of cells in the total
 * grid, number of cells and starting point of the sub grid defined
 * for one processor, cell spacing in space, time step size, etc. 
 *
 * Information about the boundary conditions for the grids is also
 * held here. This allows the subdomaining algorithm to intellegently
 * reassign boundary conditions where required. 
 *
 * This object is more convenient for a parser to work with than the
 * full grid, and it makes the subdomain algorithm easier too. 
 */
class GridInfo
{
public:
  // Global grid size (i.e. all domains). These are *sizes*; the
  // maximum index into the {e,h}{x,y,z}_ arrays is these minus one. 
  unsigned int global_dimx_;
  unsigned int global_dimy_;
  unsigned int global_dimz_;

  // Local grid starting point in global grid
  unsigned int start_x_;
  unsigned int start_y_;
  unsigned int start_z_;

  // Local grid size (this sub domain only) INCLUDING ghost cells
  // This ghost cells allow us to calculate the interior of
  // the subdomain without having to stop and ask the other
  // processors for data. 
  unsigned int dimx_;
  unsigned int dimy_;
  unsigned int dimz_;

  // These dimensions are intended to be returned by the Grid to
  // clients which need them, such as results. These sizes DO NOT
  // include the ghost cells.
  unsigned int dimx_no_sd_;
  unsigned int dimy_no_sd_;
  unsigned int dimz_no_sd_;

  // Local grid starting point W.R.T global grid, NOT INCLUDING
  // ghost cells.
  unsigned int start_x_no_sd_;
  unsigned int start_y_no_sd_;
  unsigned int start_z_no_sd_;

  // Time and space steppings; the distance between each point in the
  // grid.
  delta_t deltax_;
  delta_t deltay_;
  delta_t deltaz_;
  delta_t deltat_;

  // The default material inside the grid. Set by Grid::set_define(false)
  mat_idx_t default_mat_;

  // A grid is a cube with six faces. Those faces either need to have
  // boundary conditions, or they are subdomain boundaries and they
  // need to be shared with other processors. These arrays tell what
  // to do with each face. 
  //
  // 0 - Front (x = dimx, YZ plane)
  // 1 - Back (x = 0, YZ plane)
  // 2 - Left (y = 0, XZ plane)
  // 3 - Right (y = dimy, XZ plane)
  // 4 - Bottom (z = 0, XY plane)
  // 5 - Top (z = dimz, XY plane)
  //
protected:
  shared_ptr<BoundaryCond> face_bc_[6]; /*< Boundary condition to apply */
  
  /**
   * This region defines the size of the entire computational domain,
   * starting at zero.
   */ 
  region_t domain_;

  /**
   * This region defines the size of the computational domain present
   * on this process only, starting at zero. 
   */
  region_t local_domain_;

  /**
   * This region defines the location of the computational domain
   * running on this processor within the total computational
   * domain. The maximums will be the same as for local_domain_, but
   * the starting point will be different.
   */
  region_t local_domain_in_global_;

  /**
   * This region defines the minimums and maximums for the
   * computational domain on the local process, not including any
   * ghost cells that may exist due to the division of the global
   * computational domain among processors.
   */ 
  region_t local_domain_no_ol_;

  /**
   * This region defines the minimums and maximums for the
   * computational domain on the local process with respect to the
   * global domain, not including any ghost cells that may exist due to
   * the division of the global computational domain among processors.
   */ 
  region_t local_domain_in_global_no_ol_;

  /**
   * The order in which to apply the boundary conditions. 
   */ 
  Face bc_order_[6];

  /**
   * Compute the contents of the bc_order_ array based on the
   * following rules:
   *
   * 1) Apply E/H walls
   * 2) Apply Periodic boundaries
   * 3) Apply UPML/PML
   * 4) Apply Subdomains last
   */ 
  void reorder_boundaries();

public:
  GridInfo();

  /**
   * Copy constructor, to properly handle the dynamically allocated
   * boundary conditions. 
   *
   * @param info the GridInfo object to be copied. 
   */
  GridInfo(const GridInfo &info);

  ~GridInfo();

  /**
   * Assignment operator. To handle the dynamically allocated
   * boundary conditions properly.
   */
  GridInfo& operator=(const GridInfo &info);

  /**
   * Set the boundary condition on one of the faces of this grid. This
   * one claims ownership of the object and will delete it. 
   *
   * @param face the face to assign the boundary to. One of FRONT, BACK,
   * LEFT, RIGHT, BOTTOM, TOP as defined in Types.hh
   *
   * @param bc the boundary condition to apply. 
   *
   * @return a BoundaryCond object of the type required, in which the
   * specifics of the boundary condition can be stored.
   */ 
  void set_boundary(Face face, BoundaryCond *bc);

  /**
   * Sets the boundary condition to apply to a certian face. 
   */ 
  void set_boundary(Face face, shared_ptr<BoundaryCond> bc);  

  /**
   * Returns the type of boundary assigned to a face.
   *
   * @return BoundaryCondition from Types.hh
   */
  inline BoundaryCondition get_bc_type(Face face) const
  {
    return face_bc_[face].get()->get_type();
  }

  /** 
   * Returns a reference to the boundary condition object for a face. 
   *
   * @return ref to a BoundaryCond
   */
  inline BoundaryCond& get_boundary(Face face) const
  {
    return *(face_bc_[face].get());
  }

  /**
   * Returns the face thickness for a boundary condition
   *
   * @return an unsigned int, the thickness of the boundary condition.
   */
  unsigned int get_face_thickness(Face face) const;

  /**
   * Apply the boundary conditions to the grid. Hold the Subdomain
   * boundaries until after all the other boundaries have been
   * computed.
   *
   * @param the grid to apply to 
   */
  void apply_boundaries(Grid &grid, FieldType type);

  /**
   * Writes some cool info about the boundary conditions to the given
   * stream.
   */ 
  ostream& to_string(ostream &os) const;

};

#endif // GRID_INFO_H
