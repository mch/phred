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

#include "SurfaceCurrentResult.hh"
#include <math.h>

SurfaceCurrentResult::SurfaceCurrentResult()
  : freq_start_(0), freq_end_(0), num_freqs_(0), freq_space_(0), 
    freqs_(0), sin_(0), cos_(0)
{
  region_.xmin = region_.ymin = region_.zmin = 0;
  region_.xmax = region_.ymax = region_.zmax = 0;

  variables_["back_Jy"] = &back_Jy_;
  variables_["back_Jz"] = &back_Jz_;
  variables_["back_My"] = &back_My_;
  variables_["back_Mz"] = &back_Mz_;

  variables_["front_Jy_"] = &front_Jy_;
  variables_["front_Jz_"] = &front_Jz_;
  variables_["front_My_"] = &front_My_;
  variables_["front_Mz_"] = &front_Mz_;

  variables_["left_Jx_"] = &left_Jx_;
  variables_["left_Jz_"] = &left_Jz_;
  variables_["left_Mx_"] = &left_Mx_;
  variables_["left_Mz_"] = &left_Mz_;

  variables_["right_Jx_"] = &right_Jx_;
  variables_["right_Jz_"] = &right_Jz_;
  variables_["right_Mx_"] = &right_Mx_;
  variables_["right_Mz_"] = &right_Mz_;

  variables_["bottom_Jx_"] = &bottom_Jx_;
  variables_["bottom_Jy_"] = &bottom_Jy_;
  variables_["bottom_Mx_"] = &bottom_Mx_;
  variables_["bottom_My_"] = &bottom_My_;

  variables_["top_Jx_"] = &top_Jx_;
  variables_["top_Jy_"] = &top_Jy_;
  variables_["top_Mx_"] = &top_Mx_;
  variables_["top_My_"] = &top_My_;

}

SurfaceCurrentResult::~SurfaceCurrentResult()
{}

void init(const Grid &grid)
{
  if (box_.get())
  {
    region_ = grid.get_local_region(*(box_.get()));
  }
  
  if (do_dft_ && num_freqs_ > 2 && freq_end > freq_start_)
  {
    freq_space_ = (freq_end_ - freq_start_) / (num_freqs_ - 1);

    freqs_ = new field_t[num_freqs_];

    if (!freqs_)
      throw MemoryException();

    for (unsigned int i = 0; i < num_freqs_; i++)
    {
      freqs_ = freq_start_ + i * freq_space_;
    }
  }
  else 
    do_dft_ = false;
  
  back_Jy_.set_name(base_name_ + "_back_Jy_");
  back_Jz_.set_name(base_name_ + "_back_Jz_");
  back_My_.set_name(base_name_ + "_back_My_");
  back_Mz_.set_name(base_name_ + "_back_Mz_");

  front_Jy_.set_name(base_name_ + "_front_Jy_");
  front_Jz_.set_name(base_name_ + "_front_Jz_");
  front_My_.set_name(base_name_ + "_front_My_");
  front_Mz_.set_name(base_name_ + "_front_Mz_");

  left_Jx_.set_name(base_name_ + "_left_Jx_");
  left_Jz_.set_name(base_name_ + "_left_Jz_");
  left_Mx_.set_name(base_name_ + "_left_Mx_");
  left_Mz_.set_name(base_name_ + "_left_Mz_");

  right_Jx_.set_name(base_name_ + "_right_Jx_");
  right_Jz_.set_name(base_name_ + "_right_Jz_");
  right_Mx_.set_name(base_name_ + "_right_Mx_");
  right_Mz_.set_name(base_name_ + "_right_Mz_");

  bottom_Jx_.set_name(base_name_ + "_bottom_Jx_");
  bottom_Jy_.set_name(base_name_ + "_bottom_Jy_");
  bottom_Mx_.set_name(base_name_ + "_bottom_Mx_");
  bottom_My_.set_name(base_name_ + "_bottom_My_");

  top_Jx_.set_name(base_name_ + "_top_Jx_");
  top_Jy_.set_name(base_name_ + "_top_Jy_");
  top_Mx_.set_name(base_name_ + "_top_Mx_");
  top_My_.set_name(base_name_ + "_top_My_");

}

void deinit()
{

}

map<string, Variable *> &
SurfaceCurrentResult::get_result(const Grid &grid, 
                                 unsigned int time_step)
{
  
}
