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

#include "PlaneResult.hh"

PlaneResult::PlaneResult()
  : face_(FRONT), field_(FC_EY)
{
  plane_.x = 0;
  plane_.y = 0;
  plane_.z = 0;

  dim_lens_.push_back(0);
  dim_lens_.push_back(0);
}

PlaneResult::~PlaneResult()
{}

Data &PlaneResult::get_result(const Grid &grid, unsigned int time_step)
{
  if (result_time(time_step))
  {
    data_.set_num(1);
    data_.set_ptr(grid.get_face_start(face_, field_, plane_));

    MPI_Datatype t = grid.get_plane_dt(face_);
    data_.set_datatype(t);
  } else {
    data_.set_num(0);
  }

  return data_;
}

void PlaneResult::init(const Grid &grid)
{
  dim_lens_.clear();

  switch (face_) 
  {
  case FRONT:
  case BACK:
    dim_lens_.push_back(grid.get_ldy());
    dim_lens_.push_back(grid.get_ldz());
    break;

  case TOP:
  case BOTTOM:
    dim_lens_.push_back(grid.get_ldx());
    dim_lens_.push_back(grid.get_ldy());
    break;

  case LEFT:
  case RIGHT:
    dim_lens_.push_back(grid.get_ldx());
    dim_lens_.push_back(grid.get_ldz());
    break;
  }
}
