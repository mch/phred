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

#ifndef CSG_PRIMITIVE_H
#define CSG_PRIMITIVE_H

#include "CSGObject.hh"

/**
 * An abstract base class that represents a CSG Primitive. This just
 * provides code for dealing the the centre coordinates of all
 * primitives.
 */ 
class CSGPrimitive : public CSGObject {
public:
  CSGPrimitive();
  ~CSGPrimitive();

  virtual CSGStatus is_point_inside(float x, float y, float z) const;

  /**
   * Set the centre coordinates of the object. 
   */ 
  void set_centre(float x, float y, float z);

  /**
   * Create a copy of this object. 
   */ 
  virtual shared_ptr<CSGObject> copy() const
  {
    return shared_ptr<CSGObject>(new CSGPrimitive(*this));
  }

  /**
   * Returns the centre of the object. 
   */ 
  point get_centre() const;

  /**
   * Print a string representation to an ostream.
   */
  virtual std::ostream& to_string(std::ostream &os) const;

protected:
  float centre_[3];
};

#endif // CSG_PRIMITIVE_H
