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

#ifndef PLANE_RESULT_H
#define PLANE_RESULT_H

#include "Result.hh"
#include "../Types.hh"

#include <mpi.h>

/** 
 * A simple class that outputs the value of the field components at a
 * specific plane in space.
 */ 
class PlaneResult : public Result
{
private:
protected:
  // Position of the plane in global space. Have to translate it to
  // the local grid.
  point_t plane_;

  // The face the plane is parallel to
  Face face_;

  // The field component we are interested in. 
  FieldComponent field_;

  Variable var_; /**< Our variable */

  MPI_Datatype datatype_; 

  bool have_data_; /**< True if the node has data to contribute. */ 

public:

  PlaneResult();
  ~PlaneResult();

  /**
   * Look for the plane we want to return, and return it. 
   *
   * @param grid a reference to a Grid object
   * @return a reference to a data object, which contains an MPI
   * derived data type, a planeer, and the number of items in the
   * result.
   */
  virtual map<string, Variable *> &get_result(const Grid &grid, 
                                              unsigned int time_step);

  /**
   * Set the plane position in global coordinates
   *
   * @param p
   */
  inline void set_plane(point_t p, Face face)
  {
    plane_ = p;
    face_ = face;
  }

  /**
   * Get the position of the plane in global coordinates
   *
   * @return point_t
   */
  inline point_t get_plane()
  {
    return plane_;
  }

  /**
   * Returns the face the plane is perpendicular to. 
   * @return Face
   */
  inline Face get_face()
  {
    return face_;
  }

  /** 
   * Set the field component to return. PlaneResult returns Ey by default. 
   */
  inline void set_field(FieldComponent field)
  {
    field_ = field;
  }

  /**
   * Initalize the result. Set's the size of the plane from the grid.
   */
  virtual void init(const Grid &grid);

  /**
   * Deinit the result; free the MPI derived type. 
   */ 
  virtual void deinit();
  
  /**
   * Set the name of our variable to something human readable!
   */
  inline void set_name(const char *name)
  {
    var_.set_name(name);
  }
};

#endif // PLANE_RESULT_H
