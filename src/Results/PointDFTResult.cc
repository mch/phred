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

#include "PointDFTResult.hh"
#include "Constants.hh"
#include "Exceptions.hh"

#include <math.h>

PointDFTResult::PointDFTResult()
  : result_(0)
{
  variables_["point_DFT"] = &var_;
}

PointDFTResult::PointDFTResult(field_t freq_start,
                               field_t freq_stop, 
                               unsigned int num_freqs)
  : DFTResult(freq_start, freq_stop, num_freqs), result_(0)
{
  variables_["point_DFT"] = &var_;
}

PointDFTResult::~PointDFTResult()
{ 
  if (result_)
    delete[] result_;
}

void PointDFTResult::init(const Grid &grid)
{
  ours_ = true; 

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

  if (freq_stop_ < freq_start_)
  {
    field_t temp = freq_stop_;
    freq_stop_ = freq_start_;
    freq_start_ = temp;
  }

  var_.has_time_dimension(false); // We have only one output at the end. 
  
  var_.add_dimension("results", 13, 13, 0);
  var_.add_dimension("freqs", num_freqs_, num_freqs_, 0);
  var_.set_name(base_name_);

  freq_space_ = (freq_stop_ - freq_start_) / num_freqs_;
  result_ = new field_t[(num_freqs_ + 1) * 13];

  if (!result_)
    throw MemoryException();

  for (unsigned int i = 0; i <= num_freqs_; i++)
  {
    result_[i * 13] = freq_start_ + freq_space_ * i;
    for (unsigned int j = 1; j < 13; j++)
      result_[i*13 + j] = 0;
  }
  
  MPI_Datatype temp;
  MPI_Type_contiguous(13, GRID_MPI_TYPE, &temp);
  MPI_Type_commit(&temp);
  var_.set_datatype(temp);

  if (ours_)
  {
    var_.set_num(num_freqs_);
    var_.set_ptr(result_);
  } else {
    var_.set_num(0);
    var_.set_ptr(0);  
  }
}

void PointDFTResult::deinit(const Grid &grid)
{
  if (result_)
  {
    delete[] result_;
    result_ = 0;
  }
}

map<string, Variable *> &PointDFTResult::get_result(const Grid &grid, 
                                                     unsigned int time_step)
{
  delta_t dt = grid.get_deltat();
  delta_t time = dt * time_step;

  field_t cos_temp, sin_temp;

  if (ours_)
  {
    for (unsigned int i = 0; i <= num_freqs_; i++)
    {
      cos_temp = cos(2 * PI * result_[i*13] * time);
      sin_temp = sin(2 * PI * result_[i*13] * time);

      result_[i*13 + 1] += grid.get_ex(l_.x, l_.y, l_.z) * cos_temp;
    
      result_[i*13 + 2] += (-1) * grid.get_ex(l_.x, l_.y, l_.z) * sin_temp;

      result_[i*13 + 3] += grid.get_ey(l_.x, l_.y, l_.z) * cos_temp;
    
      result_[i*13 + 4] += (-1) * grid.get_ey(l_.x, l_.y, l_.z) * sin_temp;

      result_[i*13 + 5] += grid.get_ez(l_.x, l_.y, l_.z) * cos_temp;
    
      result_[i*13 + 6] += (-1) * grid.get_ez(l_.x, l_.y, l_.z) * sin_temp;

      // H components
      result_[i*13 + 7] += grid.get_hx(l_.x, l_.y, l_.z) * cos_temp;
    
      result_[i*13 + 8] += (-1) * grid.get_hx(l_.x, l_.y, l_.z) * sin_temp;

      result_[i*13 + 9] += grid.get_hy(l_.x, l_.y, l_.z) * cos_temp;
    
      result_[i*13 + 10] += (-1) * grid.get_hy(l_.x, l_.y, l_.z) * sin_temp;

      result_[i*13 + 11] += grid.get_hz(l_.x, l_.y, l_.z) * cos_temp;
    
      result_[i*13 + 12] += (-1) * grid.get_hz(l_.x, l_.y, l_.z) * sin_temp;
    }
  }

  return variables_;
}
