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

#include "BartlettExcitation.hh"
#include <math.h>

field_t BartlettExcitation::window(float x, float y, float z) 
{
  field_t wx = 1, wy = 1, wz = 1;

  if (xmax_ - xmin_ > 0)
    wx = 1 - fabs( ( (x - xmin_) - 0.5 * (xmax_ - xmin_)) 
                   / (0.5 * (xmax_ - xmin_)));
  
  if (ymax_ - ymin_ > 0)
    wy = 1 - fabs( ( (y - ymin_) - 0.5 * (ymax_ - ymin_)) 
                   / (0.5 * (ymax_ - ymin_)));
  
  if (zmax_ - zmin_ > 0)
    wz = 1 - fabs( ( (z - zmin_) - 0.5 * (zmax_ - zmin_)) 
                   / (0.5 * (zmax_ - zmin_)));
  
  return wx * wy * wz;
}

BartlettExcitation::BartlettExcitation(shared_ptr<SourceFunction> sf) 
  : WindowedExcitation(sf)
{}

BartlettExcitation::~BartlettExcitation()
{}
