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

#include "PointResult.hh"

PointResult::PointResult()
{
  MPI_Datatype temp;
  MPI_Type_contiguous(7, GRID_MPI_TYPE, &temp);
  MPI_Type_commit(&temp);
  data_.set_datatype(temp);

  dim_lens_.push_back(7);
}

PointResult::PointResult(point_t p)
  : point_(p)
{
  MPI_Datatype temp;
  MPI_Type_contiguous(7, GRID_MPI_TYPE, &temp);
  MPI_Type_commit(&temp);
  data_.set_datatype(temp);

  dim_lens_.push_back(7);
}

PointResult::~PointResult()
{ }

Data &PointResult::get_result(const Grid &grid, unsigned int time_step)
{
  point_t l;
  bool ours = true;

  if (result_time(time_step))
  {
    if (point_.x < grid.get_lsx() 
        || point_.x >= grid.get_lsx() + grid.get_ldx())
      ours = false;
    else 
      l.x = point_.x - grid.get_lsx();

    if (point_.y < grid.get_lsy() 
        || point_.y >= grid.get_lsy() + grid.get_ldy())
      ours = false;
    else 
      l.y = point_.y - grid.get_lsy();

    if (point_.z < grid.get_lsz() 
        || point_.z >= grid.get_lsz() + grid.get_ldz())
      ours = false;
    else 
      l.z = point_.z - grid.get_lsz();

    if (ours) 
    {
      field_data_[0] = grid.get_deltat() * time_step;
      field_data_[1] = grid.get_ex(l.x, l.y, l.z);
      field_data_[2] = grid.get_ey(l.x, l.y, l.z);
      field_data_[3] = grid.get_ez(l.x, l.y, l.z);
    
      field_data_[4] = grid.get_hx(l.x, l.y, l.z);
      field_data_[5] = grid.get_hy(l.x, l.y, l.z);
      field_data_[6] = grid.get_hz(l.x, l.y, l.z);
    
      data_.set_ptr(field_data_);
      data_.set_num(1);
    } 
    else
    {
      data_.set_ptr(0);
      data_.set_num(0);
    }
  } else {
    data_.set_num(0);
  }

  return data_;
}
