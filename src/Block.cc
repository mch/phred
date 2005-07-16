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

#include "Block.hh"

std::ostream &operator<<(std::ostream &os, const Block &b)
{
  if (b.is_global())
  {
    os << "Block of grid cells, in the global grid.";
  } else {
    os << "Block of grid cells, in the local grid.";
  }
  
  os << "\n\tLower lefthand "
    "corner is at (" << b.xmin() << ", " << b.ymin() << ", " << b.zmin()
            << "), extends to\n\tand includes ("
            << b.xmax() << ", " << b.ymax() << ", " << b.zmax() 
            << ").\n\tLength: " << b.xlen() << ", " << b.ylen() 
            << ", " << b.zlen()
            << ".";

  if (!b.is_global())
  {
    os << "\n\tWith respect to the block of grid cells in the "
      "global domain,\n\tthis local block starts at, or is offset by, (" 
       << b.xoffset() << ", " << b.yoffset() << ", " << b.zoffset()
       << "). ";
  }

  if (b.has_data())
  {
    os << "\n\tThis block contains data. " << std::endl;

//     os << "This block has faces [";
    
//     os << " in the grid. " << endl;
  } else {
    os << "\n\tThis block DOES NOT contain data. " << std::endl;
  }

  return os;
}
