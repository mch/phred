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

    for (unsigned int i = 0; i < num_freqs_; i++)
    {
      freqs_ = freq_start_ + i * freq_space_;
    }
  }
  else 
    do_dft_ = false;
  
  
}

void deinit()
{

}

map<string, Variable *> &
SurfaceCurrentResult::get_result(const Grid &grid, 
                                 unsigned int time_step)
{

}
