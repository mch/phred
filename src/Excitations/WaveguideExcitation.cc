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

#include "WaveguideExcitation.hh"
#include "../Constants.hh"
#include <math.h>

WaveguideExcitation::WaveguideExcitation(shared_ptr<SourceFunction> sf)
  : WindowedExcitation(sf)
{}

WaveguideExcitation::~WaveguideExcitation()
{}

field_t WaveguideExcitation::window(float x, float y, float z)
{
  unsigned int dx, dy, dz;
  field_t ret = 1.0;

  if (box_.get())
  {
    point size = (*box_).get_size();
  
    if (size.x > 0 && mode_x_ > 0)
      ret = ret * sin(mode_x_ * PI * (x - xmin_) / size.x);
    
    if (size.y > 0 && mode_y_ > 0)
      ret = ret * sin(mode_y_ * PI * (y - ymin_) / size.y);
    
    if (size.z > 0 && mode_z_ > 0)
      ret = ret * sin(mode_z_ * PI * (z - zmin_) / size.z);
    
    // if (y == 10)
    //     cerr << "Waveguide excitation at (" << x << ", " << y << ", " 
    //          << z << ") is " << ret
    //          << endl;
  }

  return ret;
}

ostream& WaveguideExcitation::to_string(ostream &os) const
{
  os << "WaveguideExcitation with modes " << mode_x_
     << ", " << mode_y_ << ", and " << mode_z_ 
     << " being applied to a region starting at ("
     << xmin_ << ", " << ymin_ << ", " << zmin_ 
     << ") and extending to (" << xmax_ << ymax_
     << zmax_ << ". On this rank, the excitation starts at ("
     << lxmin_ << ", " << lymin_ << ", " << lzmin_ 
     << "). ";

  Excitation::to_string(os);

  return os;    
}
