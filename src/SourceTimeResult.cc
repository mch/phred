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

#include "SourceTimeResult.hh"

SourceTimeResult::SourceTimeResult(SourceFunction &te)
  : te_(te)
{
  dim_lens_.push_back(2);
  var_name_ = "Source Time Excitation";

  MPI_Datatype temp;
  MPI_Type_contiguous(2, GRID_MPI_TYPE, &temp);
  MPI_Type_commit(&temp);

  data_.set_ptr(result_);
  data_.set_datatype(temp);
}

SourceTimeResult::~SourceTimeResult()
{ }

void SourceTimeResult::set_excitation(const SourceFunction &te)
{
  te_ = te;
}

Data &SourceTimeResult::get_result(const Grid &grid, unsigned int time_step)
{
  if (result_time(time_step)) 
  {
    result_[0] = grid.get_deltat() * time_step; 
    result_[1] = te_.source_function(grid, time_step);
    data_.set_num(1); 
  } 
  else {
    data_.set_num(0); 
  }

  return data_;
}
