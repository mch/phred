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

#include "Gaussm.hh"

Gaussm::Gaussm()
  : alpha_(1), deltaf_(0.5), f0_(1)
{}

Gaussm::~Gaussm()
{}

field_t Gaussm::signal_function(float time) const
{
  field_t temp = (time - 4. / (PI * deltaf_)) * deltaf_ * PI;

  return alpha_ * exp((-1) * temp * temp)
    * sin(2. * PI * f0_ * (time - 4. / (PI * deltaf_)));
}

void Gaussm::set_parameters(field_t alpha, field_t deltaf, field_t f0)
{
  alpha_ = alpha;
  deltaf_ = deltaf;
  f0_ = f0;
}

field_t Gaussm::get_alpha() const
{
  return alpha_;
}

field_t Gaussm::get_deltaf() const
{
  return deltaf_;
}

field_t Gaussm::get_f0() const
{
  return f0_;
}

field_t Gaussm::length() const
{
  return 8.0 / (PI * deltaf_);
}

ostream& Gaussm::to_string(ostream &os) const
{
  return os << "Gaussm signal with an amplitude of " << alpha_
            << ", a centre frequency of " << f0_
            << " Hz, and a width of " << deltaf_ << " Hz. Length: "
            << length() << " seconds.";
}
