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

#include "CSGCylinder.hh"

CSGCylinder::CSGCylinder()
  : radius_(0.5), height_(1)
{}

CSGCylinder::~CSGCylinder()
{}

CSGStatus CSGCylinder::is_point_inside(float x, float y, float z) const
{
  CSGStatus ret = OUTSIDE;

  float r = sqrt(pow(centre_[0] - x, 2) + pow(centre_[1] - y, 2));

  if (r < radius_ && z > centre_[2] - height_ / 2 
      && z < centre_[2] + height_ / 2)
    ret = INSIDE;

  return ret; 
}

void CSGCylinder::set_radius(float radius)
{
  if (radius > 0)
    radius_ = radius;
}

float CSGCylinder::get_radius() const
{
  return radius_;
}

void CSGCylinder::set_height(float height)
{
  if (height > 0)
    height_ = height;
}
  
float CSGCylinder::get_height() const
{
  return height_;
}
