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

#include "BlockResult.hh"

BlockResult::BlockResult()
  : field_comp_(FC_EY), init_(false)
{
  
}

BlockResult::BlockResult(region_t r, FieldComponent field_comp)
  : region_(r), field_comp_(field_comp), init_(false)
{
}

BlockResult::~BlockResult()
{

}

void BlockResult::init(const Grid &grid)
{
  MPI_Datatype temp;
  int sizes[3];
  int subsizes[3];
  int starts[3];

  // Setup (convert to local)
  region_ = grid.global_to_local(region_);
  sizes[0] = grid.get_ldx();
  sizes[1] = grid.get_ldy();
  sizes[2] = grid.get_ldz();

  subsizes[0] = region_.xmax - region_.xmin;
  subsizes[1] = region_.ymax - region_.ymin;
  subsizes[2] = region_.zmax - region_.zmin;

  starts[0] = region_.xmin;
  starts[1] = region_.ymin;
  starts[2] = region_.zmin;

  // Create
  MPI_Type_create_subarray(3, sizes, subsizes, starts, 1, 
                           GRID_MPI_TYPE, &temp);
  MPI_Type_commit(&temp);
  data_.set_datatype(temp);
  data_.set_num(0);
  data_.set_ptr(const_cast<field_t *>(grid.get_pointer(point_t(region_.xmin, 
                                                               region_.ymin, 
                                                               region_.zmin), 
                                                       field_comp_)));
  init_ = true; 
}

Data &BlockResult::get_result(const Grid &grid, unsigned int time_step)
{
  if (init_ && result_time(time_step))
  {
    data_.set_num(1);
  } 
  else 
  {
    data_.set_num(0);
  }

  return data_;
}
