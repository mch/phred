/* 
   Phred - Phred is a parallel finite difference time domain
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

#include "GaussPulse.hh"

GaussPulse::GaussPulse()
  : alpha_(1), tau_p_(1)
{}

GaussPulse::~GaussPulse()
{}

field_t GaussPulse::signal_function(float time) const
{
  field_t temp = time / tau_p_ - 4;

  return alpha_ * exp(-1 * temp * temp / 2);
}

void GaussPulse::set_parameters(field_t alpha, field_t tau_p)
{
  alpha_ = alpha;
  tau_p_ = tau_p;
}

field_t GaussPulse::get_alpha() const
{
  return alpha_;
}

field_t GaussPulse::get_taup() const
{
  return tau_p_;
}

field_t GaussPulse::length() const
{
  return 8 * tau_p_;
}

ostream& GaussPulse::to_string(ostream &os) const
{
  return os << "GaussPulse signal with an amplitude of " << alpha_
            << ", a characteristic time of " << tau_p_
            << ", and a length of "
            << length() << " seconds.";
}
