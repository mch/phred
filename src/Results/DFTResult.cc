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

#include "DFTResult.hh"

DFTResult::DFTResult()
  : freq_start_(0), freq_stop_(0), num_freqs_(0)
{}

DFTResult::DFTResult(field_t freq_start, field_t freq_stop, 
                     unsigned int num_freqs)
  : freq_start_(freq_start), freq_stop_(freq_stop), num_freqs_(num_freqs)
{}

DFTResult::~DFTResult()
{}

void DFTResult::set_freq(field_t freq_start, field_t freq_stop, 
                         unsigned int num_freqs)
{
  if (freq_start > freq_stop)
    throw ResultException("Starting frequency must be less than stop frequency for DFTResults.");

  if (num_freqs < 2)
    throw ResultException("DFTResults must return at least two frequencies.");

  freq_start_ = freq_start;
  freq_stop_ = freq_stop;
  num_freqs_ = num_freqs;
}
