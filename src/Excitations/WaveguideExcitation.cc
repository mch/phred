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

#include "WaveguideExcitation.hh"
#include "Constants.hh"
#include <math.h>

WaveguideExcitation::WaveguideExcitation(SourceFunction *sf)
  : WindowedExcitation(sf)
{}

WaveguideExcitation::~WaveguideExcitation()
{}

field_t WaveguideExcitation::window(region_t r, 
                                    unsigned int x, 
                                    unsigned int y, 
                                    unsigned int z)
{
  unsigned int dx, dy, dz;
  field_t ret = 1.0;

  dx = region_.xmax - region_.xmin;
  dy = region_.ymax - region_.ymin;
  dz = region_.zmax - region_.zmin;

  if (dx > 0 && mode_x_ > 0)
    ret = ret * sin(mode_x_ * PI * (x - region_.xmin) / dx);

  if (dy > 0 && mode_y_ > 0)
    ret = ret * sin(mode_y_ * PI * (y - region_.ymin) / dy);
  
  if (dz > 0 && mode_z_ > 0)
    ret = ret * sin(mode_z_ * PI * (z - region_.zmin) / dz);
  
  // if (y == 10)
//     cerr << "Waveguide excitation at (" << x << ", " << y << ", " 
//          << z << ") is " << ret
//          << endl;

  return ret;
}

