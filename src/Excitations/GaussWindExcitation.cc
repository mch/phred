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

GaussWindExcitation::GaussWindExcitation(SourceFunction *sf)
  : WindowedExcitation(sf)
{}

GaussWindExcitation::~GaussWindExcitation()
{}

field_t GaussWindExcitation::window(float x, float y, float z) 
{
  field_t wx = 1, wy = 1, wz = 1;

  point centre = bx_.get_centre();

  float x = (i - r.xmin) - 0.5 * (r.xmax - r.xmin - 1); 
  float y = (j - r.ymin) - 0.5 * (r.ymax - r.ymin - 1); 
  float z = (k - r.zmin) - 0.5 * (r.zmax - r.zmin - 1); 

  wx = exp((-pow(x, 2))/(2*pow(std_dev_, 2)));
  wy = exp((-pow(y, 2))/(2*pow(std_dev_, 2)));
  wz = exp((-pow(z, 2))/(2*pow(std_dev_, 2)));

  return wx * wy * wz;
}

void GaussWindExcitation::set_std_dev(float stddev)
{
  std_dev_ = stddev;
}

float GaussWindExcitation::get_std_dev()
{
  return std_dev_;
}
