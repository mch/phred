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

#ifndef CSG_CYLINDER_H
#define CSG_CYLINDER_H

#include "CSGPrimitive.hh"

/**
 * A cylinder, which is oriented along the z axis, having a default
 * height of 1 and a default radius of 0.5.
 */ 
class CSGCylinder : public CSGPrimitive 
{
public:
  CSGCylinder();
  ~CSGCylinder();

  /**
   * Tell us if a point is inside or outside the solid.
   *
   * @param x x coordinate of the point
   * @param y y coordinate of the point
   * @param z z coordinate of the point
   * @return 1 if the point is inside the solid, 0 if the point is on
   * the boundary, and -1 if it is outside. 
   */
  CSGStatus is_point_inside(float x, float y, float z) const;

  /**
   * Create a copy of this object. 
   */ 
  virtual shared_ptr<CSGObject> copy() const
  {
    return shared_ptr<CSGObject>(new CSGCylinder(*this));
  }

  /**
   * Set the radius of this cylinder. 
   */ 
  void set_radius(float radius);

  /**
   * Returns the radius of this cylinder
   */ 
  float get_radius() const;

  /**
   * Set the height of this cylinder.
   */
  void set_height(float height);
  
  /**
   * Returns the height of this cylinder
   */ 
  float get_height() const;

  /**
   * Print a string representation to an ostream.
   */
  virtual std::ostream& to_string(std::ostream &os) const;

protected:
  float radius_;
  float height_;

};

#endif // CSG_CYLINDER_H
