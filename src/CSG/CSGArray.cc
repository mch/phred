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
  : xspace_(0), yspace_(0), zspace_(0), xlen_(1), ylen_(1), zlen_(1), 
    child_(obj)
{}

CSGArray::CSGArray(const CSGArray &rhs)
{
  *this = rhs;
}

const CSGArray &CSGArray::operator= (const CSGArray &rhs)
{
  child_ = (*rhs.child_).copy();

  return *this;
}

CSGArray::~CSGArray()
{}

CSGStatus CSGArray::is_point_inside(float x, float y, float z) const
{
  CSGStatus ret = OUTSIDE; 
  float tx, ty, tz;

  tx = x; 
  // This tricky part! I think this can be done better:
  for (unsigned int i = 0; i < xlen_; i++, tx = tx - xspace_)
  {
    ty = y;
    for (unsigned int j = 0; j < ylen_; j++, ty = ty - yspace_)
    {
      tz = z;
      for (unsigned int k = 0; k < zlen_; k++, tz = tz - zspace_)
      {
        ret = (*child_).is_point_inside(tx, ty, tz);
        if (ret != OUTSIDE)
          return ret;
      }
    }
  }

  return ret;
}

shared_ptr<CSGObject> CSGArray::copy() const
{
  return shared_ptr<CSGObject> (new CSGArray(*this));
}

void CSGArray::set_lengths(unsigned int xlen, unsigned int ylen, 
                           unsigned int zlen)
{
  if (xlen < 1)
    throw CSGException("The length of the array along the x dimension must be greater than or equal to 1.");

  if (ylen < 1)
    throw CSGException("The length of the array along the y dimension must be greater than or equal to 1.");

  if (zlen < 1)
    throw CSGException("The length of the array along the z dimension must be greater than or equal to 1.");

  xlen_ = xlen;
  ylen_ = ylen;
  zlen_ = zlen;
}

void CSGArray::set_spacing(float x, float y, float z)
{
  if (x < 0)
    throw CSGException("x spacing must not be less than zero.");

  if (y < 0)
    throw CSGException("y spacing must not be less than zero.");

  if (z < 0)
    throw CSGException("z spacing must not be less than zero.");

  xspace_ = x;
  yspace_ = y;
  zspace_ = z;
}

std::ostream& CSGArray::to_string(std::ostream &os) const
{
  return os << "CSGArray, array is " << xlen_ << " x " << ylen_ << " x "
            << xlen_ << ", with array factors of " << xspace_ 
            << " by " << yspace_ << " by " << zspace_;
}
