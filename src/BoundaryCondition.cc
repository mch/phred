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

#include "BoundaryCondition.hh"
#include "Grid.hh" // Flesh out the forward declaration

// PML's are actually one cell thicker than advertised, and this
// routine bears that out...
region_t BoundaryCond::find_face(Face face, Grid &grid)
{
  region_t r;

  switch (face) 
  {
  case FRONT:
    r.xmax = grid.get_ldx();
    r.xmin = r.xmax - 1 - thickness_;
    //r.xmin = r.xmax - thickness_;
    r.ymin = r.zmin = 0;
    r.ymax = grid.get_ldy();
    r.zmax = grid.get_ldz();
    break;

  case BACK:
    r.xmax = 1 + thickness_;
    //r.xmax = thickness_;
    r.xmin = 0;
    r.ymin = r.zmin = 0;
    r.ymax = grid.get_ldy();
    r.zmax = grid.get_ldz();
    break;

  case TOP:
    r.zmin = grid.get_ldz() - 1 - thickness_;
    //r.zmin = grid.get_ldz() - thickness_;
    r.zmax = grid.get_ldz();
    r.ymin = r.xmin = 0;
    r.xmax = grid.get_ldx();
    r.ymax = grid.get_ldy();
    break;

  case BOTTOM:
    r.zmin = 0;
    r.zmax = 1 + thickness_;
    //r.zmax = thickness_;
    r.ymin = r.xmin = 0;
    r.xmax = grid.get_ldx();
    r.ymax = grid.get_ldy();
    break;

  case LEFT:
    r.ymin = 0;
    r.ymax = 1 + thickness_;
    //r.ymax = thickness_;
    r.zmin = r.xmin = 0;
    r.xmax = grid.get_ldx();
    r.zmax = grid.get_ldz();
    break;

  case RIGHT:
    r.ymin = grid.get_ldy() - 1 - thickness_;
    //r.ymin = grid.get_ldy() - thickness_;
    r.ymax = grid.get_ldy();
    r.zmin = r.xmin = 0;
    r.xmax = grid.get_ldx();
    r.zmax = grid.get_ldz();
    break;
  }

  return r; 
}

// void BoundaryCond::apply(Face face, Grid &grid)
// {
//   region_t r = find_face(face, grid);

//   switch (face)
//   {
//   case FRONT:
//   case BACK:
//     condition<YZPlane>(r, grid);
//     break;

//   case LEFT:
//   case RIGHT:
//     condition<XZPlane>(r, grid);
//     break;

//   case TOP:
//   case BOTTOM:
//     helper<XYPlane>(r, grid);
//     break;
//   }
  
// }

unsigned int BoundaryCond::get_thickness()
{
  return thickness_;
}
