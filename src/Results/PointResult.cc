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
#include "../Globals.hh"
#include "../Types.hh"

PointResult::PointResult()
{
  variables_["Point result"] = &var_;
}

PointResult::PointResult(grid_point p)
  : point_(p)
{
  variables_["Point result"] = &var_;
}

PointResult::~PointResult()
{ }

void PointResult::calculate_result(const Grid &grid, 
                                   unsigned int time_step)
{
  if (ours_ && result_time(time_step)) 
  {
    field_data_[0] = grid.get_deltat() * time_step;
    field_data_[1] = grid.get_ex(l_.x, l_.y, l_.z);
    field_data_[2] = grid.get_ey(l_.x, l_.y, l_.z);
    field_data_[3] = grid.get_ez(l_.x, l_.y, l_.z);
    
    field_data_[4] = grid.get_hx(l_.x, l_.y, l_.z);
    field_data_[5] = grid.get_hy(l_.x, l_.y, l_.z);
    field_data_[6] = grid.get_hz(l_.x, l_.y, l_.z);
    
    var_.set_num(1);
  } 
  else
  {
    var_.set_num(0);
  }
  
}

// MCH, 2005-02-09: changed > to >= for lines like 
// point_.y >= grid.get_lsy() + grid.get_ldy()) because
// grid.get_lsy() + grid.get_ldy()) is one over the maximum indexable
// cell in the local region.
void PointResult::init(const Grid &grid)
{
  ours_ = true; 

  point_ = grid.get_global_cell(space_point_);

  if (point_.x < grid.get_lsx() 
      || point_.x >= grid.get_lsx() + grid.get_ldx())
    ours_ = false;
  else 
    l_.x = point_.x - grid.get_lsx();

  if (point_.y < grid.get_lsy() 
      || point_.y >= grid.get_lsy() + grid.get_ldy())
    ours_ = false;
  else 
    l_.y = point_.y - grid.get_lsy();

  if (point_.z < grid.get_lsz() 
      || point_.z >= grid.get_lsz() + grid.get_ldz())
    ours_ = false;
  else 
    l_.z = point_.z - grid.get_lsz();

#ifdef DEBUG
  if (ours_)
    cerr << "PointResult at " << l_.x << "x" << l_.y
         << "x" << l_.z << " belongs to " << MPI_RANK << endl;
  else
    cerr << "PointResult at " << l_.x << "x" << l_.y
         << "x" << l_.z << " DOES NOT belong to " << MPI_RANK << endl;
#endif

  MPI_Datatype temp;
  MPI_Type_contiguous(7, GRID_MPI_TYPE, &temp);
  MPI_Type_commit(&temp);
  var_.set_datatype(temp);

  var_.add_dimension("field components", 7, 7, 0);
  var_.set_name(base_name_);
  var_.set_ptr(field_data_);
}

void PointResult::deinit()
{}

ostream& PointResult::to_string(ostream &os) const
{
  os << "PointResult: sampling data at point " << space_point_
     << "\nsampling data at grid cell " << point_ << "\n";

  return os;
}
