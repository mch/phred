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

#include "FakeResult.hh"
#include "../Globals.hh"

FakeResult::FakeResult() 
  : data_(0)
{}
FakeResult::~FakeResult() 
{
  if (data_)
    delete[] data_;
}

void FakeResult::init(const Grid &grid)
{
  data_ = new field_t[2 * 8];
  
  for (int i = 0; i < 16; i++)
    data_[i] = MPI_RANK + i;

  MPI_Type_contiguous(16, GRID_MPI_TYPE, &dtype_);
  MPI_Type_commit(&dtype_);

  var_.add_dimension("x", 2, 2 * MPI_SIZE, 2 * MPI_RANK);
  var_.add_dimension("y", 8, 8, 0);

  var_.has_time_dimension(false);

  variables_["fake"] = &var_;
  var_.set_name(base_name_);
  var_.set_ptr(data_);
  var_.set_datatype(dtype_);
}

void FakeResult::deinit()
{
  MPI_Type_free(&dtype_);

  if (data_)
  {
    delete[] data_;
    data_ = 0;
  }
}

map<string, Variable *> &
FakeResult::get_result(const Grid &grid, 
                       unsigned int time_step)
{
  if (result_time(time_step))
  {
    var_.set_num(1);
  } else {
    var_.set_num(0);
  }
  
  return variables_;
}

ostream& FakeResult::to_string(ostream &os) const
{
  return os << "Fake result for DataWriter testing.";
}
