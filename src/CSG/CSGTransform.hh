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

#ifndef CSG_TRANSFORM_H
#define CSG_TRANSFORM_H

#include "CGSObject.hh"

/**
 * This class transforms a child object by rotation, translation,
 * and scaling.
 */ 
class CSGTransform {
public:
  CSGTransform(const CGSObject &child);
  virtual ~CSGTransform();
  
  /**
   * Tell us if a point is inside or outside the solid.
   *
   * @param x x coordinate of the point
   * @param y y coordinate of the point
   * @param z z coordinate of the point
   * @return true if the point is inside the solid.
   */
  bool is_point_inside(float x, float y, float z) const;
  
  /**
   * Returns true if the object on the right hand side is enclosed by
   * the transform on the left.
   */
  bool operator> (const CSGTransform &rhs) const;

  /**
   * Returns true if the object on the left hand side is enclosed by the 
   * transform on the left.
   */ 
  bool operator< (const CSGTransform &rhs) const;

  /**
   * Returns a list of verticies for this object. The
   */ 
  const vector<float[3]> &get_verticies() const;

protected:
  const CSGObject &child_;
};

#endif // CSG_TRANSFORM_H
