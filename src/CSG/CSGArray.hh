/* 
   Phred - Phred is a parallel finite difference time domain
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

#ifndef CSG_ARRAY_H
#define CSG_ARRAY_H

#include "CSGObject.hh"

/**
 * This class replicates a child object into an array along at least
 * one dimension.
 */ 
class CSGArray : public CSGObject
{
public:
  CSGArray(shared_ptr<CSGObject> obj);

  CSGArray(const CSGArray &rhs);

  const CSGArray &operator= (const CSGArray &rhs);

  ~CSGArray();

  /**
   * Tell us if a point is inside or outside the solid.
   *
   * @param x x coordinate of the point
   * @param y y coordinate of the point
   * @param z z coordinate of the point
   * @return the status of the point, inside, outside, or on the surface.
   */
  CSGStatus is_point_inside(float x, float y, float z) const;
  
  /**
   * Create a copy of this object. 
   */ 
  shared_ptr<CSGObject> copy() const;

private:
  // Spacing between the elements in each direction
  float xspace_, yspace_, zspace_;

  // Length of the array in each direction. Must be greater than one. 
  unsigned int lenx_, leny_, lenz_;

  // The child object we are making into an array. 
  shared_ptr<CSGObject> child_;
};

#endif
