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

#include "Block.hh"

ostream &Block::operator<<(ostream &os)
{
  return os << "Block of grid cells, in the local grid, lower lefthand "
    "corner is at (" << xmin_ << ", " << ymin_ << ", " << zmin_
            << "), extends to but does not include ("
            << xmax_ ", " << ymax_ << ", " << zmax_ << "). Length: "
            << len_x_ << ", " << len_y_ << ", " << len_z_ 
            << ". With respect to the block of grid cells in the "
    "global domain, this local block starts at (" 
            << start_x_ << ", " << start_y_ << ", " << start_z_
            << "). ";
}
