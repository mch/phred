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

#ifndef CSG_BOX_H
#define CSG_BOX_H

#include "CSGPrimitive.hh"

/**
 * A box in three dimensions, centred at 0,0,0 and 1,1,1 in size by
 * default. 
 */ 
class CSGBox : public CSGPrimitive 
{
public:
  CSGBox();
  ~CSGBox();

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
   * Set the size of the box by specifying the lengths along the x, y
   * and z axis.
   */ 
  void set_size(float x_size, float y_size, float z_size);

  /**
   * Returns the size of the box
   */ 
  point get_size() const;

  /**
   * Create a copy of this object. 
   */ 
  virtual shared_ptr<CSGObject> copy() const
  {
    return shared_ptr<CSGObject>(new CSGBox(*this));
  }

  /**
   * Returns true if the object on the right hand side is enclosed by
   * the object on the left.
   */
  //bool operator> (const CSGObject &rhs) const;

  /**
   * Returns true if the object on the left hand side is enclosed by the 
   * object on the left.
   */ 
  //bool operator< (const CSGObject &rhs) const;

  /**
   * Returns a list of verticies for this object.
   */ 
  //const vector<float[3]> &get_verticies() const;

  /**
   * Print a string representation to an ostream.
   */
  virtual std::ostream& to_string(std::ostream &os) const;

private:
  float lengths_[3];
  float half_lengths_[3];
};

#endif // CSG_BOX_H
