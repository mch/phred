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

#include "SphereGeom.hh"
#include "Grid.hh"

#include <math.h>

Sphere::Sphere()
{}

Sphere::~Sphere()
{}

void Sphere::init(const Grid &grid)
{

}

void Sphere::set_material(Grid &grid)
{
  region_t r;

  int temp = centre_.x - radius_;
  r.xmin = temp >= 0 ? temp : 0;
  temp = centre_.x + radius_;
  r.xmax = temp < grid.get_ldx_sd() 
    ? temp : 0;

  temp = centre_.y - radius_;
  r.ymin = (centre_.y - radius_) >= 0 ? (centre_.y - radius_) : 0;
  temp = centre_.y + radius_;
  r.ymax = temp < grid.get_ldy_sd() ? temp : 0;

  temp = centre_.z - radius_;
  r.zmin = temp >= 0 ? temp : 0;
  temp = centre_.z + radius_;
  r.zmax = temp < grid.get_ldz_sd() ? temp : 0;

  for (unsigned int i = r.xmin; i < r.xmax; i++)
  {
    for (unsigned int j = r.ymin; j < r.ymax; j++)
    {
      // oops, fix this:
      int p = radius_ * radius_ - (i - centre_.x) * (i - centre_.x) 
        - (j - centre_.y) * (j - centre_.y);

      if (p > 0)
      {
        unsigned int zmax = 
          static_cast<unsigned int>(centre_.z + 
                                    ceil(sqrt(static_cast<double>(p))));
        unsigned int zmin = 
          static_cast<unsigned int>(centre_.z - 
                                    ceil(sqrt(static_cast<double>(p))));

        for (unsigned int k = zmin; k < zmax; k++)
        {
          grid.set_material(i, j, k, material_id_);
        }
      }
    }
  }
}
