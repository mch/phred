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

#include "Region.hh"

Region::Region()
  : global_xmin_(0), global_xmax_(0),
    global_ymin_(0), global_ymax_(0),
    global_zmin_(0), global_zmax_(0),
    global_(0)
{}

Region::Region(const Region &global)
  : global_xmin_(0), global_xmax_(0),
    global_ymin_(0), global_ymax_(0),
    global_zmin_(0), global_zmax_(0)
{
  global_ = new Region(global);
}

Region::Region(unsigned int xmin, unsigned int xmax,
               unsigned int ymin, unsigned int ymax,
               unsigned int zmin, unsigned int zmax)
  : global_xmin_(xmin), global_xmax_(xmax),
    global_ymin_(ymin), global_ymax_(ymax),
    global_zmin_(zmin), global_zmax_(zmax)
{}

Region::Region(unsigned int xmin, unsigned int xmax,
               unsigned int ymin, unsigned int ymax,
               unsigned int zmin, unsigned int zmax,
               const Region &global)
{
  global_ = new Region(global);
}

Region::~Region()
{
  if (global_)
    delete global_;
}

bool Region::is_local() const
{
  bool ret = false;
  if (global_xmin_ >= grid_info.start_x_no_sd_ 
      && global_xmax_ < grid_info.start_x_no_sd_ + grid_info.dimx_no_sd_
      && global_ymin_ >= grid_info.start_y_no_sd_ 
      && global_ymax_ < grid_info.start_y_no_sd_ + grid_info.dimy_no_sd_
      && global_zmin_ >= grid_info.start_z_no_sd_ 
      && global_zmax_ < grid_info.start_z_no_sd_ + grid_info.dimz_no_sd_)
    ret = true;

  return ret;
}

void Region::set_x(unsigned int xmin, unsigned int xmax)
{
  if (xmin > xmax)
    throw RegionException("xmin cannot be greater than xmax!");

//   if (xmax > grid_info.global_dimx_)
//     throw RegionException("xmax is outside of the grid's region.");

  global_xmin_ = xmin;
  global_xmax_ = xmax;
}

void Region::set_y(unsigned int ymin, unsigned int ymax)
{
  if (ymin > ymax)
    throw RegionException("ymin cannot be greater than ymax!");

//   if (ymax > grid_info.global_dimy_)
//     throw RegionException("ymax is outside of the grid's region.");

  global_ymin_ = ymin;
  global_zmax_ = ymax;
}

void Region::set_z(unsigned int zmin, unsigned int zmax)
{
  if (zmin > zmax)
    throw RegionException("zmin cannot be greater than zmax!");

//   if (zmax > grid_info.global_dimz_)
//     throw RegionException("zmax is outside of the grid's region.");

  global_zmin_ = zmin;
  global_zmax_ = zmax;
}

unsigned int Region::get_xmin() const
{
  unsigned int ret = 0;
  if (global_xmin_ < grid_info.start_x_no_sd_
      && global_xmax_ >= grid_info.start_x_no_sd_)
    ret = grid_info.start_x_no_sd_;
  else if (global_xmin_ >= grid_info.start_x_no_sd_
           && global_xmin_ < grid_info.start_x_no_sd_ 
                             + grid_info.dimx_no_sd_)
    ret = global_xmin_ - grid_info.start_x_no_sd_;

  return ret;
}

unsigned int Region::get_ymin() const
{
  unsigned int ret = 0;
  if (global_ymin_ < grid_info.start_y_no_sd_
      && global_ymax_ >= grid_info.start_y_no_sd_)
    ret = grid_info.start_y_no_sd_;
  else if (global_ymin_ >= grid_info.start_y_no_sd_
           && global_ymin_ < grid_info.start_y_no_sd_ 
                             + grid_info.dimy_no_sd_)
    ret = global_ymin_ - grid_info.start_y_no_sd_;

  return ret;
}

unsigned int Region::get_zmin() const
{
  unsigned int ret = 0;
  if (global_zmin_ < grid_info.start_z_no_sd_
      && global_zmax_ >= grid_info.start_z_no_sd_)
    ret = grid_info.start_z_no_sd_;
  else if (global_zmin_ >= grid_info.start_z_no_sd_
           && global_zmin_ < grid_info.start_z_no_sd_ 
                             + grid_info.dimz_no_sd_)
    ret = global_zmin_ - grid_info.start_z_no_sd_;

  return ret;
}

unsigned int Region::get_xmax() const
{
  unsigned int upper_bnd = grid_info.start_x_no_sd_ 
    + grid_info.dimx_no_sd_ - 1;
  unsigned int ret = 0;

  if (global_xmax_ >= upper_bnd && global_xmin_ <= upper_bnd)
    ret = grid_info.dimx_no_sd_ - 1;
  else if (global_xmax_ >= grid_info.start_x_no_sd_)
    ret = global_xmax_ - grid_info.start_x_no_sd_;

  return ret;
}

unsigned int Region::get_ymax() const
{
  unsigned int upper_bnd = grid_info.start_y_no_sd_ 
    + grid_info.dimy_no_sd_ - 1;
  unsigned int ret = 0;

  if (global_ymax_ >= upper_bnd && global_ymin_ <= upper_bnd)
    ret = grid_info.dimy_no_sd_ - 1;
  else if (global_ymax_ >= grid_info.start_y_no_sd_)
    ret = global_ymax_ - grid_info.start_y_no_sd_;

  return ret;
}

unsigned int Region::get_zmax() const
{
  unsigned int upper_bnd = grid_info.start_z_no_sd_ 
    + grid_info.dimz_no_sd_ - 1;
  unsigned int ret = 0;

  if (global_zmax_ >= upper_bnd && global_zmin_ <= upper_bnd)
    ret = grid_info.dimz_no_sd_ - 1;
  else if (global_zmax_ >= grid_info.start_z_no_sd_)
    ret = global_zmax_ - grid_info.start_z_no_sd_;

  return ret;
}

///////////////////////////////////////////////////////////////
// OverlapRegion implementation
///////////////////////////////////////////////////////////////

OverlapRegion::OverlapRegion(const Region &region)
  : Region(region)
{}

OverlapRegion::OverlapRegion(unsigned int xmin, unsigned int xmax,
                             unsigned int ymin, unsigned int ymax,
                             unsigned int zmin, unsigned int zmax, 
                             const Region &region)
  : Region(xmin, xmax, ymin, ymax, zmin, zmax, region)
{}

OverlapRegion::~OverlapRegion()
{}

bool OverlapRegion::is_local() const
{
  bool ret = false;
  if (global_xmin_ >= grid_info.start_x_ 
      && global_xmax_ < grid_info.start_x_ + grid_info.dimx_
      && global_ymin_ >= grid_info.start_y_ 
      && global_ymax_ < grid_info.start_y_ + grid_info.dimy_
      && global_zmin_ >= grid_info.start_z_ 
      && global_zmax_ < grid_info.start_z_ + grid_info.dimz_)
    ret = true;

  return ret;
}

unsigned int OverlapRegion::get_xmin() const
{
  unsigned int ret = 0;
  if (global_xmin_ < grid_info.start_x_
      && global_xmax_ >= grid_info.start_x_)
    ret = grid_info.start_x_;
  else if (global_xmin_ >= grid_info.start_x_
           && global_xmin_ < grid_info.start_x_ 
                             + grid_info.dimx_)
    ret = global_xmin_ - grid_info.start_x_;

  return ret;
}

unsigned int OverlapRegion::get_ymin() const
{
  unsigned int ret = 0;
  if (global_ymin_ < grid_info.start_y_
      && global_ymax_ >= grid_info.start_y_)
    ret = grid_info.start_y_;
  else if (global_ymin_ >= grid_info.start_y_
           && global_ymin_ < grid_info.start_y_ 
                             + grid_info.dimy_)
    ret = global_ymin_ - grid_info.start_y_;

  return ret;
}

unsigned int OverlapRegion::get_zmin() const
{
  unsigned int ret = 0;
  if (global_zmin_ < grid_info.start_z_
      && global_zmax_ >= grid_info.start_z_)
    ret = grid_info.start_z_;
  else if (global_zmin_ >= grid_info.start_z_
           && global_zmin_ < grid_info.start_z_ 
                             + grid_info.dimz_)
    ret = global_zmin_ - grid_info.start_z_;

  return ret;
}

unsigned int OverlapRegion::get_xmax() const
{
  unsigned int upper_bnd = grid_info.start_x_ 
    + grid_info.dimx_ - 1;
  unsigned int ret = 0;

  if (global_xmax_ >= upper_bnd && global_xmin_ <= upper_bnd)
    ret = grid_info.dimx_ - 1;
  else if (global_xmax_ >= grid_info.start_x_)
    ret = global_xmax_ - grid_info.start_x_;

  return ret;
}

unsigned int OverlapRegion::get_ymax() const
{
  unsigned int upper_bnd = grid_info.start_y_ 
    + grid_info.dimy_ - 1;
  unsigned int ret = 0;

  if (global_ymax_ >= upper_bnd && global_ymin_ <= upper_bnd)
    ret = grid_info.dimy_ - 1;
  else if (global_ymax_ >= grid_info.start_y_)
    ret = global_ymax_ - grid_info.start_y_;

  return ret;
}

unsigned int OverlapRegion::get_zmax() const
{
  unsigned int upper_bnd = grid_info.start_z_ 
    + grid_info.dimz_ - 1;
  unsigned int ret = 0;

  if (global_zmax_ >= upper_bnd && global_zmin_ <= upper_bnd)
    ret = grid_info.dimz_ - 1;
  else if (global_zmax_ >= grid_info.start_z_)
    ret = global_zmax_ - grid_info.start_z_;

  return ret;
}
