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

#include "PeriodicExcitation.hh"

PeriodicExcitation::PeriodicExcitation(shared_ptr<Signal> sf)
  : Excitation(sf)
{}

PeriodicExcitation::~PeriodicExcitation()
{}

void PeriodicExcitation::excite(Grid &grid, 
                                unsigned int time_step, 
                                FieldType type)
{
  Excitation::excite(grid, time_step, type);
}

void PeriodicExcitation::init(const Grid &grid)
{
  
  
}

ostream& PeriodicExcitation::to_string(ostream &os) const
{
  os << "PeriodicExcitation... " << endl;
  
  return os;
}

void PeriodicExcitation::set_region(shared_ptr<CSGBox> box, Face face)
{
  face_ = face;

  // Make sure that the box is actually 0 cells thick on the face we
  // want to apply it, make a copy if necessary.
  
  Excitation::set_region(box);
}

void PeriodicExcitation::set_poynting(float x, float y, float z)
{
  // GNDN for now. Once this is enabled, the excite() function will
  // need a lot of work!
}
