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

#ifndef SPHERE_GEOM_H_
#define SPHERE_GEOM_H_

#include "Geometry.hh"

/**
 * Defines a sphere in the grid
 */
class Sphere : public Geometry
{
private:
protected:
  point_t centre_;
  unsigned int radius_;

public:
  Sphere();
  Sphere(point_t centre, unsigned int radius);
  ~Sphere();

  /**
   * Set the radius of the sphere. 
   */
  inline void set_radius(unsigned int r)
  {
    radius_ = r;
  }

  /**
   * Set the centre of the sphere. 
   */
  inline void set_centre(point_t c)
  {
    centre_ = c;
  }

  /**
   * Init the region, set the material in the grid. Normally called by
   * the FDTD object.
   */
  virtual void init(const Grid &grid);
  
  /**
   * Update the material indicies of the grid 
   */
  virtual void set_material(Grid &grid);

  /**
   * Returns true if the given point is inside this geometry, with
   * respect to local grid coordinates. 
   */ 
  virtual bool local_point_inside(unsigned int x,
                                  unsigned int y, 
                                  unsigned int z);
};

#endif
