/* 
   Phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004-2005 Matt Hughes <mhughe@uvic.ca>

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
#include <cmath>

using namespace std;

CSGCylinder::CSGCylinder()
  : radius_(0.5), height_(1)
{}

CSGCylinder::~CSGCylinder()
{}

CSGStatus CSGCylinder::is_point_inside(float x, float y, float z) const
{
  CSGStatus ret = OUTSIDE;

  if (z > centre_[2] - height_ / 2 
      && z < centre_[2] + height_ / 2
      && y > centre_[1] - radius_
      && y < centre_[1] + radius_
      && x > centre_[0] - radius_
      && x < centre_[0] + radius_)
  {
    float temp1 = centre_[0] - x;
    float temp2 = centre_[1] - y;
    //float r = sqrt(temp1*temp1 + temp2*temp2);
    float r_sq = temp1*temp1 + temp2*temp2;

    //if (r < radius_)
    if (r_sq < (radius_ * radius_))
      ret = INSIDE;
  }

  return ret; 
}

void CSGCylinder::set_radius(float radius)
{
  if (radius > 0)
    radius_ = radius;
  else
    throw CSGException("Cylinder radius must be greater than zero.");
}

float CSGCylinder::get_radius() const
{
  return radius_;
}

void CSGCylinder::set_height(float height)
{
  if (height > 0)
    height_ = height;
  else
    throw CSGException("Cylinder height must be greater than zero.");
}
  
float CSGCylinder::get_height() const
{
  return height_;
}

std::ostream& CSGCylinder::to_string(std::ostream &os) const
{
  return os << "CSGCylinder, centred at (" << centre_[0] << ", " << centre_[1]
            << ", " << centre_[2] << "), with a radius of " << radius_
            << " and a height of " << height_;
}
