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

#include "CSGIntersection.hh"

CSGIntersection::CSGIntersection(shared_ptr<CSGObject> left, 
                                 shared_ptr<CSGObject> right)
  : CSGOperator(left, right)
{}

CSGIntersection::~CSGIntersection()
{}

CSGStatus CSGIntersection::is_point_inside(float x, float y, float z) const
{
  CSGStatus r = (*right_).is_point_inside(x, y, z);
  CSGStatus l = (*left_).is_point_inside(x, y, z);

  CSGStatus ret = l;

  if (l == INSIDE && r == INSIDE)
    ret = INSIDE;

  return ret;
}

std::ostream& CSGIntersection::to_string(std::ostream &os) const
{
  if (right_.get() && left_.get())
  {
    os << "CSGIntersection, taking the intersection of a " 
       << (*right_) << " and " << " a " << (*left_);
  } else {
    os << "CSGIntersection, missing objects";
  }
  return os;
}
