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

#include "CSGBox.hh"

CSGBox::CSGBox()
{
  set_size(1, 1, 1);
}

CSGBox::~CSGBox()
{}

CSGStatus CSGBox::is_point_inside(float x, float y, float z) const
{
  CSGStatus ret = OUTSIDE;

  if (x > centre_[0] - half_lengths_[0] 
      && x < centre_[0] + half_lengths_[0]
      && y > centre_[1] - half_lengths_[1] 
      && y < centre_[1] + half_lengths_[1]
      && z > centre_[2] - half_lengths_[2] 
      && z < centre_[2] + half_lengths_[2])
    ret = INSIDE;
  
  if ((x == centre_[0] - half_lengths_[0] 
       || x == centre_[0] + half_lengths_[0])
      && (y == centre_[1] - half_lengths_[1] 
          || y == centre_[1] + half_lengths_[1])
      && (z == centre_[2] - half_lengths_[2] 
          || z < centre_[2] + half_lengths_[2]))
    ret = BOUNDARY;

  return ret;
}

void CSGBox::set_size(float x_size, float y_size, float z_size)
{
  lengths_[0] = x_size;
  lengths_[1] = y_size;
  lengths_[2] = z_size;

  half_lengths_[0] = x_size / 2.0;
  half_lengths_[1] = y_size / 2.0;
  half_lengths_[2] = z_size / 2.0;

}

point CSGBox::get_size()
{
  return point(lengths_[0], lengths_[1], lengths_[2]);
}
