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

#include "SourceDFTResult.hh"
#include "Constants.hh"

#include <math.h>

SourceDFTResult::SourceDFTResult(SourceFunction &te)
  : te_(te), freq_start_(0), freq_stop_(0),
    num_freqs_(0)
{}

SourceDFTResult::SourceDFTResult(SourceFunction &te, 
                                  field_t freq_start,
                                  field_t freq_stop, 
                                  unsigned int num_freqs)
  : te_(te), freq_start_(freq_start), freq_stop_(freq_stop),
    num_freqs_(num_freqs)
{}

SourceDFTResult::~SourceDFTResult()
{ 
  if (result_)
    delete[] result_;
}

void SourceDFTResult::init(const Grid &grid)
{
  if (freq_stop_ < freq_start_)
  {
    field_t temp = freq_stop_;
    freq_stop_ = freq_start_;
    freq_start_ = temp;
  }

  time_dim_ = false; // We have only one output at the end. 

  dim_lens_.push_back(3); // Rows of (Freq, DFT)
  dim_lens_.push_back(num_freqs_); // num_freq rows
  var_name_ = "Source DFT";

  freq_space_ = (freq_stop_ - freq_start_) / num_freqs_;
  result_ = new field_t[(num_freqs_ + 1) * 3];

  for (unsigned int i = 0; i <= num_freqs_; i++)
  {
    result_[i * 3] = freq_start_ + freq_space_ * i;
    result_[i*3 + 1] = 0;
    result_[i*3 + 2] = 0;
  }
  
  MPI_Datatype temp;
  MPI_Type_contiguous(3, GRID_MPI_TYPE, &temp);
  MPI_Type_commit(&temp);

  data_.set_num(num_freqs_);
  data_.set_ptr(result_);
  data_.set_datatype(temp);
}

void SourceDFTResult::deinit(const Grid &grid)
{
  if (result_)
    delete[] result_;
}

void SourceDFTResult::set_excitation(const SourceFunction &te)
{
  te_ = te;
}

Data &SourceDFTResult::get_result(const Grid &grid, unsigned int time_step)
{
  field_t sf = te_.source_function(grid, time_step);
  for (unsigned int i = 0; i <= num_freqs_; i++)
  {
    result_[i*3 + 1] += sf * cos(2 * PI * result_[i*3] 
                                 * time_step * grid.get_deltat());

    result_[i*3 + 2] += (-1) * sf * sin(2 * PI * result_[i*3] 
                                        * time_step * grid.get_deltat());
  }

  return data_;
}
