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

#ifndef CSG_OBJECT_H
#define CSG_OBJECT_H

#include "../Types.hh"
#include "../Exceptions.hh"

#include <boost/shared_ptr.hpp>

using namespace boost;

/**
 * This is an abstract base class for a node in a constructive 
 * solid geometry binary tree. 
 */ 
class CSGObject {
public:
  virtual ~CSGObject()
  {}

  /**
   * Tell us if a point is inside or outside the solid.
   *
   * @param x x coordinate of the point
   * @param y y coordinate of the point
   * @param z z coordinate of the point
   * @return the status of the point, inside, outside, or on the surface.
   */
  virtual CSGStatus is_point_inside(float x, float y, float z) const
  { return OUTSIDE; } // Should be abstract; this makes boost.python happy. 
  
  /**
   * Create a copy of this object. 
   */ 
  virtual shared_ptr<CSGObject> copy() const
  {
    return shared_ptr<CSGObject> (new CSGObject(*this));
  }

  /**
   * Returns true if the object on the right hand side is enclosed by
   * the object on the left.
   */
  //virtual bool operator> (const CSGObject &rhs) const = 0;

  /**
   * Returns true if the object on the left hand side is enclosed by the 
   * object on the left.
   */ 
  //virtual bool operator< (const CSGObject &rhs) const = 0;

  /**
   * Returns a list of verticies for this object. The
   */ 
  //virtual const vector<float[3]> &get_verticies() const = 0;


  /**
   * Print a string representation to an std::ostream.
   */
  virtual std::ostream& to_string(std::ostream &os) const
  {
    return os << "A CSGObject of indeterminate type.";
  }

  friend std::ostream& operator<< (std::ostream& os, const CSGObject &c);
};

inline std::ostream& operator<< (std::ostream& os, const CSGObject &c)
{
  return c.to_string(os);
}

#endif // CSG_OBJECT_H
