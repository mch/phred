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

#include "CSGArray.hh"

CSGArray::CSGArray(shared_ptr<CSGObject> obj)
  : xspace_(0), yspace_(0), zspace_(0), lenx_(1), leny_(1), lenz_(1), 
    child_(obj)
{}

CSGArray::CSGArray(const CSGArray &rhs)
{
  *this = rhs;
}

const CSGArray &CSGArray::operator= (const CSGArray &rhs)
{
  child_ = (*rhs.child).copy();

  return *this;
}

CSGArray::~CSGArray()
{}

CSGStatus CSGArray::is_point_inside(float x, float y, float z) const
{
  // This tricky part!
}

shared_ptr<CSGObject> CSGArray::copy() const
{
  return shared_ptr<CSGObject> (new CSGArray(*this));
}
