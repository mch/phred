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

#include "FarfieldResult.hh"

FarfieldResult::FarfieldResult()
  : theta_start_(0), theta_stop_(0), phi_start_(90), phi_stop_(90),
    axis_(X_AXIS), freq_start_(0), freq_stop_(100), num_freqs_(10), 
    freq_space_(0), e_theta_re_(0), e_theta_im_(0), e_phi_re_(0), 
    e_phi_im_(0), jff_mom_(0), mff_mom_(0)
{}

FarfieldResult~FarfieldResult()
{}

Data & FarfieldResult::get_result(const Grid &grid, unsigned int time_step)
{
  
}

void FarfieldResult::init(const Grid &grid)
{

}
  
void FarfieldResult::deinit(const Grid &grid)
{

}

