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

BoundaryCond::BoundaryCond() 
  : thickness_(0) 
{}

BoundaryCond::~BoundaryCond() 
{}

// PML's are actually one cell thicker than advertised, and this
// routine bears that out...
region_t BoundaryCond::find_face(Face face, const Grid &grid)
{
  region_t r;

  switch (face) 
  {
  case FRONT:
    r.xmax = grid.get_ldx_sd();
    r.xmin = r.xmax - 1 - thickness_;
    //r.xmin = r.xmax - thickness_;
    r.ymin = r.zmin = 0;
    r.ymax = grid.get_ldy_sd();
    r.zmax = grid.get_ldz_sd();
    break;

  case BACK:
    r.xmax = 1 + thickness_;
    //r.xmax = thickness_;
    r.xmin = 0;
    r.ymin = r.zmin = 0;
    r.ymax = grid.get_ldy_sd();
    r.zmax = grid.get_ldz_sd();
    break;

  case TOP:
    r.zmin = grid.get_ldz_sd() - 1 - thickness_;
    //r.zmin = grid.get_ldz() - thickness_;
    r.zmax = grid.get_ldz_sd();
    r.ymin = r.xmin = 0;
    r.xmax = grid.get_ldx_sd();
    r.ymax = grid.get_ldy_sd();
    break;

  case BOTTOM:
    r.zmin = 0;
    r.zmax = 1 + thickness_;
    //r.zmax = thickness_;
    r.ymin = r.xmin = 0;
    r.xmax = grid.get_ldx_sd();
    r.ymax = grid.get_ldy_sd();
    break;

  case LEFT:
    r.ymin = 0;
    r.ymax = 1 + thickness_;
    //r.ymax = thickness_;
    r.zmin = r.xmin = 0;
    r.xmax = grid.get_ldx_sd();
    r.zmax = grid.get_ldz_sd();
    break;

  case RIGHT:
    r.ymin = grid.get_ldy_sd() - 1 - thickness_;
    //r.ymin = grid.get_ldy() - thickness_;
    r.ymax = grid.get_ldy_sd();
    r.zmin = r.xmin = 0;
    r.xmax = grid.get_ldx_sd();
    r.zmax = grid.get_ldz_sd();
    break;
  }

  return r; 
}

void BoundaryCond::compute_regions(Face face, const Grid &grid)
{
  grid_r_ = find_face(face, grid);

  bc_r_.xmin = bc_r_.ymin = bc_r_.zmin = 0;

  bc_r_.xmax = grid_r_.xmax - grid_r_.xmin;
  bc_r_.ymax = grid_r_.ymax - grid_r_.ymin;
  bc_r_.zmax = grid_r_.zmax - grid_r_.zmin;

  grid_ex_r_ = grid_ey_r_ = grid_ez_r_ = grid_r_;
  grid_hx_r_ = grid_hy_r_ = grid_hz_r_ = grid_r_;

  grid_ex_r_.ymin++;
  grid_ex_r_.zmin++;
  
  grid_ey_r_.xmin++;
  grid_ey_r_.zmin++;
  
  grid_ez_r_.xmin++;
  grid_ez_r_.ymin++;
  
  grid_hx_r_.ymax--;
  grid_hx_r_.zmax--;
  
  grid_hy_r_.xmax--;
  grid_hy_r_.zmax--;
  
  grid_hz_r_.xmax--;
  grid_hz_r_.ymax--;

  // Corrections made by comparing to Jan's FDTD
  grid_ex_r_.xmax--;
  grid_ey_r_.ymax--;
  grid_ez_r_.zmax--;

  // Don't allow external face E field updates (electric walls)
  // Make sure that the PML computes all components at internal faces. 
  if (thickness_ > 0) 
  {
    switch (face) {
    case FRONT:
      grid_ey_r_.xmin--;
      grid_ez_r_.xmin--;

      grid_ey_r_.xmax--;
      grid_ez_r_.xmax--;

      grid_ex_r_.ymax--;
      grid_ez_r_.ymax--;

      grid_ex_r_.zmax--;
      grid_ey_r_.zmax--;
      break;

    case BACK:
      //grid_hz_r_.xmax++;
      //grid_hy_r_.xmax++;

      grid_ex_r_.ymax--;
      grid_ez_r_.ymax--;

      grid_ex_r_.zmax--;
      grid_ey_r_.zmax--;
      break;

    case TOP:
      grid_ex_r_.zmin--;
      grid_ey_r_.zmin--;

      grid_ex_r_.zmax--;
      grid_ey_r_.zmax--;

      grid_ey_r_.xmax--;
      grid_ez_r_.xmax--;

      grid_ez_r_.ymax--;
      grid_ex_r_.ymax--;
      break;

    case BOTTOM:
      //grid_hy_r_.zmax++;
      //grid_hx_r_.zmax++;

      grid_ey_r_.xmax--;
      grid_ez_r_.xmax--;

      grid_ez_r_.ymax--;
      grid_ex_r_.ymax--;
      break;

    case LEFT:
      //grid_hz_r_.ymax++;
      //grid_hx_r_.ymax++;

      grid_ey_r_.xmax--;
      grid_ez_r_.xmax--;

      grid_ex_r_.zmax--;
      grid_ey_r_.zmax--;
      break;

    case RIGHT:
      grid_ez_r_.ymin--;
      grid_ex_r_.ymin--;

      grid_ez_r_.ymax--;
      grid_ex_r_.ymax--;

      grid_ey_r_.xmax--;
      grid_ez_r_.xmax--;

      grid_ex_r_.zmax--;
      grid_ey_r_.zmax--;
      break;
    }
  }
}

region_t BoundaryCond::find_local_region(region_t field_r)
{
  region_t r;

// The right way: 
  r.xmin = field_r.xmin - grid_r_.xmin;
  r.ymin = field_r.ymin - grid_r_.ymin;
  r.zmin = field_r.zmin - grid_r_.zmin;
  
  r.xmax = r.xmin  + (field_r.xmax - field_r.xmin);
  r.ymax = r.ymin  + (field_r.ymax - field_r.ymin);
  r.zmax = r.zmin  + (field_r.zmax - field_r.zmin);
  
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

void BoundaryCond::set_thickness(unsigned int thickness)
{
  thickness_ = thickness;
}

unsigned int BoundaryCond::get_thickness() const
{
  return thickness_;
}

BoundaryCondition BoundaryCond::get_type() const
{
  return UNKNOWN;
}

void BoundaryCond::init(const Grid &grid, Face face)
{}

void BoundaryCond::deinit(const Grid &grid, Face face)
{}

void BoundaryCond::add_sd_bcs(SubdomainBc *sd, Face bcface, Face sdface)
{}
