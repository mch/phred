/* 
   phred - Phred is a parallel finite difference time domain
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

#ifndef POINT_RESULT_H
#define POINT_RESULT_H

#include "Result.hh"
#include "../Types.hh"

#include <mpi.h>

/** 
 * A simple class that outputs the value of the field components at a
 * specific point in space.
 */ 
class PointResult : public Result
{
private:
protected:
  /** 
   * The field data the we copy from the grid. In order: time, ex, ey, ez,
   * hx, hy, hz.
   */
  field_t field_data_[7];

  point space_point_; /**< Point in real space, as set by the user */ 

  grid_point point_; /**< Point in global space. Have to translate it to
                     the local grid. */

  grid_point l_; /**< Point in local space. */

  // True if the point is in our part of the grid
  bool ours_;

  // Our variable
  Variable var_;

public:

  PointResult();
  PointResult(grid_point p);
  ~PointResult();

  /**
   * Look for the point we want to return, and return it. 
   *
   * @param grid a reference to a Grid object
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  void calculate_result(const Grid &grid, 
                        unsigned int time_step);

  /**
   * Set the point in global coordinates
   *
   * @param point
   */
  inline void set_point(point p)
  {
    space_point_ = p;
  }

  /**
   * Get the point in global coordinates
   *
   * @return grid_point
   */
  inline point get_point()
  {
    return space_point_;
  }
  
  void init(const Grid &grid);
  void deinit();

  /**
   * Print a string representation to an ostream.
   */
  ostream& to_string(ostream &os) const;

};

#endif // POINT_RESULT_H
