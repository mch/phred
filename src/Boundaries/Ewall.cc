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

#include "Ewall.hh"
//#include "../GridPlane.hh"
#include "../Grid.hh"

Ewall::Ewall()
{}

Ewall::~Ewall()
{}

// template <class T>
// void Ewall::condition(region_t r, Grid &grid)
// {
//   T p(grid);

//   for (unsigned int i = r.xmin; i < r.xmax; i++) {
//     for (unsigned int j = r.ymin; j < r.ymax; j++) {
//       for (unsigned int k = r.zmin; k < r.zmax; k++) {
//         p.set_e_t1(i, j, k, 0.);
//         p.set_e_t2(i, j, k, 0.);
//       }
//     }
//   }
// }

void Ewall::apply(Face face, Grid &grid, FieldType type)
{
  if (type != E && type != BOTH)
    return;
  
  region_t r = find_face(face, grid);

  switch (face)
  {
  case FRONT:
    for (int j = 0; j < grid.get_ldy_sd(); j++)
    {
      for (int k = 0; k < grid.get_ldz_sd(); k++)
      {
        grid.set_ey(grid.get_ldx_sd() - 1, j, k, 0.0);
        grid.set_ez(grid.get_ldx_sd() - 1, j, k, 0.0);
      }
    }
    break;

  case BACK:
    for (int j = 0; j < grid.get_ldy_sd(); j++)
    {
      for (int k = 0; k < grid.get_ldz_sd(); k++)
      {
        grid.set_ey(0, j, k, -1 * grid.get_ey(1, j, k));
        grid.set_ez(0, j, k, -1 * grid.get_ez(1, j, k));
      }
    }

    //condition<YZPlane>(r, grid);
    break;

  case LEFT:
    for (int i = 0; i < grid.get_ldx_sd(); i++)
    {
      for (int k = 0; k < grid.get_ldz_sd(); k++)
      {
        grid.set_ex(i, 0, k, -1 * grid.get_ex(i, 1, k));
        grid.set_ez(i, 0, k, -1 * grid.get_ez(i, 1, k));
      }
    }
    break;

  case RIGHT:
    for (int i = 0; i < grid.get_ldx_sd(); i++)
    {
      for (int k = 0; k < grid.get_ldz_sd(); k++)
      {
        grid.set_ex(i, grid.get_ldy_sd() - 1, k, 0.0);
        grid.set_ez(i, grid.get_ldy_sd() - 1, k, 0.0);
      }
    }
    //condition<XZPlane>(r, grid);
    break;

  case TOP:
    for (int i = 0; i < grid.get_ldx_sd(); i++)
    {
      for (int j = 0; j < grid.get_ldy_sd(); j++)
      {
        grid.set_ex(i, j, grid.get_ldz_sd() - 1, 0.0);
        grid.set_ey(i, j, grid.get_ldz_sd() - 1, 0.0);
      }
    }
    break;

  case BOTTOM:
    for (int i = 0; i < grid.get_ldx_sd(); i++)
    {
      for (int j = 0; j < grid.get_ldy_sd(); j++)
      {
        grid.set_ex(i, j, 0, -1 * grid.get_ex(i, j, 1));
        grid.set_ey(i, j, 0, -1 * grid.get_ey(i, j, 1));
      }
    }
    //condition<XYPlane>(r, grid);
    break;
  }
  
}

BoundaryCondition Ewall::get_type() const
{
  return EWALL;
}
