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

#ifndef BOX_H_
#define BOX_H_

#include "Geometry.hh"

/**
 * Defines a box in the grid
 */
class Box : public Geometry
{
private:
protected:

public:
  Box();
  Box(region_t r);
  ~Box();

  /**
   * Set the region to operate on
   */
  inline void set_region(region_t r)
  {
    bounding_box_ = r;
  }

  /**
   * Set the region to operate on
   */
  void set_region(unsigned int xstart, unsigned int xstop, 
                  unsigned int ystart, unsigned int ystop, 
                  unsigned int zstart, unsigned int zstop);

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
