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

#include "CSGObject.hh"

/**
 * This class transforms a child object by rotation, translation,
 * and scaling. x_t = A * x_0 + t
 *
 * \bug Scaling is not implemented
 * \bug rotation about an arbitrary point is not implemented. 
 */ 
class CSGTransform : public CSGObject {
public:
  CSGTransform(shared_ptr<CSGObject> child);
  ~CSGTransform();
  
  /**
   * Tell us if a point is inside or outside the solid.
   *
   * @param x x coordinate of the point
   * @param y y coordinate of the point
   * @param z z coordinate of the point
   * @return true if the point is inside the solid.
   */
  CSGStatus is_point_inside(float x, float y, float z) const;
  
  /**
   * Create a copy of this object. 
   */ 
  virtual shared_ptr<CSGObject> copy() const;

  /**
   * Set up rotation
   *
   * @param The point to rotate about
   * @param The axis to rotate around
   * @param The angle to rotate
   */ 
  void set_rotation(point p, point vector, float angle);

  /**
   * Set up scaling
   */ 
  void set_scaling(float sx, float sy, float sz);

  /**
   * Set up translation
   */ 
  void set_translation(float tx, float ty, float tz);

protected:
  // The child being transformed
  shared_ptr<CSGObject> child_; 

  // Translation
  float tx_, ty_, tz_;
  
  // Rotation about a point using quaternions
  point vector_;
  float angle_;
  point rotation_point_;

  // Scaling about a point
  float sx_, sy_, sz_;

  // Precalculated transformation matrix for scaling and rotation
  float A[3][3];

  // Precalculated inverse transformation 
  float Ainv[3][3];

  /**
   * Calculate a rotation matrix. A must be a 3x3 array. 
   */ 
  void calc_rotation_matrix(const point &v, 
                            const float &angle, 
                            float A[3][3]) const;
};

#endif // CSG_TRANSFORM_H
