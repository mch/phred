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

#include "Box.hh"
#include "Grid.hh"

Box::Box()
{}

Box::~Box()
{}

void Box::set_region(unsigned int xstart, unsigned int xstop, 
                     unsigned int ystart, unsigned int ystop, 
                     unsigned int zstart, unsigned int zstop)
{
  r_.xmin = xstart;
  r_.xmax = xstop;
  r_.ymin = ystart;
  r_.ymax = ystop;
  r_.zmin = zstart;
  r_.zmax = zstop;
}

void Box::init(const Grid &grid)
{

}
void Box::set_material(Grid &grid)
{
  region_t r = grid.global_to_local(r_);

  for (unsigned int i = r.xmin; i < r.xmax; i++)
  {
    for (unsigned int j = r.ymin; j < r.ymax; j++)
    {
      for (unsigned int k = r.zmin; k < r.zmax; k++)
      {
        grid.set_material(i, j, k, material_id_);
      }
    }
  }
}
