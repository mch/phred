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

#include "GridInfo.hh"

GridInfo::GridInfo() 
  : global_dimx_(0), global_dimy_(0), global_dimz_(0), 
    start_x_(0), start_y_(0), start_z_(0),
    dimx_(0), dimy_(0), dimz_(0), 
    dimx_no_sd_(0), dimy_no_sd_(0), dimz_no_sd_(0), 
    start_x_no_sd_(0), start_y_no_sd_(0), start_z_no_sd_(0),
    deltax_(0), deltay_(0), deltaz_(0), deltat_(0)
{
  // Default the boundaries to ewalls
  for (int i = 0; i < 6; i++)
  {
    face_bc_[i] = shared_ptr<BoundaryCond>(new Ewall());
  }
}

GridInfo::GridInfo(const GridInfo &info) {
  *this = info;
}

GridInfo::~GridInfo()
{

}

GridInfo& GridInfo::operator=(const GridInfo &info)
{
  if (this == &info) return *this;

  global_dimx_ = info.global_dimx_;
  global_dimy_ = info.global_dimy_;
  global_dimz_ = info.global_dimz_;
  
  start_x_ = info.start_x_;
  start_y_ = info.start_y_;
  start_z_ = info.start_z_;
  
  dimx_ = info.dimx_;
  dimy_ = info.dimy_;
  dimz_ = info.dimz_;

  dimx_no_sd_ = info.dimx_no_sd_;
  dimy_no_sd_ = info.dimy_no_sd_;
  dimz_no_sd_ = info.dimz_no_sd_;

  start_x_no_sd_ = info.start_x_no_sd_;
  start_y_no_sd_ = info.start_y_no_sd_;
  start_z_no_sd_ = info.start_z_no_sd_;
  
  deltax_ = info.deltax_;
  deltay_ = info.deltay_;
  deltaz_ = info.deltaz_;
  deltat_ = info.deltat_;
  
  for (int i = 0; i < 6; i++) {
    face_bc_[i] = info.face_bc_[i];
  }
  
  return *this;
}

void GridInfo::set_boundary(Face face, BoundaryCond *bc)
{
  face_bc_[face] = shared_ptr<BoundaryCond>(bc);
}

void GridInfo::set_boundary(Face face, shared_ptr<BoundaryCond> bc)
{
    face_bc_[face] = bc;
}

unsigned int GridInfo::get_face_thickness(Face face) const
{
  unsigned int ret = 0;

  if (face_bc_[face].get()) {
    ret = face_bc_[face].get()->get_thickness();
  }

  return ret;
}

void GridInfo::apply_boundaries(Grid &grid, FieldType type)
{
  for (unsigned int i = 0; i < 6; i++)
  {
    if (face_bc_[i].get() && face_bc_[i].get()->get_type() != SUBDOMAIN) {
      face_bc_[i].get()->apply(static_cast<Face>(i), grid, type);
    }
  }

  for (unsigned int i = 0; i < 6; i++)
  {
    if (face_bc_[i].get() && face_bc_[i].get()->get_type() == SUBDOMAIN) {
      face_bc_[i].get()->apply(static_cast<Face>(i), grid, type);
    }
  }
}

