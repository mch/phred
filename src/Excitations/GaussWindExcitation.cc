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

#include "GaussWindExcitation.hh"
#include <math.h>

GaussWindExcitation::GaussWindExcitation(shared_ptr<SourceFunction> sf)
  : WindowedExcitation(sf)
{}

GaussWindExcitation::~GaussWindExcitation()
{}

field_t GaussWindExcitation::window(float x, float y, float z) 
{
  field_t wx = 1, wy = 1, wz = 1;

  if (box_.get())
  {
    point centre = (*box_).get_centre();

    float i = (x - xmin_) - centre.x;
    float j = (y - ymin_) - centre.y;
    float k = (z - zmin_) - centre.z;

    wx = exp((-(i*i)) / sdev_x_);
    wy = exp((-(j*j)) / sdev_y_);
    wz = exp((-(k*k)) / sdev_z_);
  }

  return wx * wy * wz;
}

void GaussWindExcitation::init(const Grid &grid)
{
  WindowedExcitation::init(grid);

  sdev_x_ = (xmax_ - xmin_) * 0.1;
  sdev_y_ = (ymax_ - ymin_) * 0.1;
  sdev_z_ = (zmax_ - zmin_) * 0.1;
}
