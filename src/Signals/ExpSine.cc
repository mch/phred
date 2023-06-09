/* 
   phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004-2005 Matt Hughes <mhughe@uvic.ca>

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

#include "ExpSine.hh"
#include "../Constants.hh"

ExpSine::ExpSine()
  : period_(1), omega_(1*2*PI), ampl_(1)
{}

ExpSine::ExpSine(float frequency)
  : period_(1/frequency), omega_(frequency * 2 * PI), ampl_(1)
{}

ExpSine::~ExpSine()
{}

field_t ExpSine::signal_function(float time) const
{
  field_t ret = ampl_ * (1 - exp(-1 * time / period_))
    * sin(omega_ * time);

  return ret;
}

ostream& ExpSine::to_string(ostream &os) const
{
  return os << "ExpSine signal with an amplitude of " << ampl_
            << ", a sine frequency of " << omega_ / (2 * PI) 
            << "Hz, and a period of " << period_ << " seconds.";
}
