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

#include "CSGSphere.hh"
#include <math.h>

CSGSphere::CSGSphere()
  : radius_(0.5)
{}

CSGSphere::~CSGSphere()
{}

CSGStatus CSGSphere::is_point_inside(float x, float y, float z) const
{
  CSGStatus ret = OUTSIDE;

  float r = sqrt(pow(centre_[0] - x, 2) + pow(centre_[1] - y, 2) 
                 + pow(centre_[2] - z, 2));

  if (r < radius_)
    ret = INSIDE;
  else if (r == radius_)
    ret = BOUNDARY;

  return ret; 
}

void CSGSphere::set_radius(float radius)
{
  radius_ = radius;
}

float CSGSphere::get_radius() const
{
  return radius_;
}
