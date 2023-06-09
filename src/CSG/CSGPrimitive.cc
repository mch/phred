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

#include "CSGPrimitive.hh"

CSGPrimitive::CSGPrimitive()
{
  centre_[0] = 0;
  centre_[1] = 0;
  centre_[2] = 0;  
}

CSGPrimitive::~CSGPrimitive()
{}

// This is just here to make boost.python happy. 
CSGStatus CSGPrimitive::is_point_inside(float x, float y, float z) const
{
  return OUTSIDE;
}

void CSGPrimitive::set_centre(float x, float y, float z)
{
  centre_[0] = x;
  centre_[1] = y;
  centre_[2] = z;
}

point CSGPrimitive::get_centre() const
{
  return point(centre_[0], centre_[1], centre_[2]);
}

std::ostream& CSGPrimitive::to_string(std::ostream &os) const
{
  return os << "A CSGPrimitive of indeterminate type";
}


