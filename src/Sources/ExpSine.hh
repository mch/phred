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

#ifndef EXCITE_EXPSINE_H
#define EXCITE_EXPSINE_H

#include "SourceFunction.hh"
#include "../Constants.hh"

#include <math.h>

class ExpSine : public SourceFunction
{
public:
  /**
   * Produces a sine function that is windows by an exponential
   * function so that the initial discontinuity is not so bad. 
   *
   * @param time the time at which to apply the excitation
   *
   * @return the value of the excitation
   */
  field_t source_function(float time);

  ExpSine();
  ExpSine(float frequency);
  ~ExpSine();

  /**
   * Set the frequency of the excitation in Hz. 
   */
  inline void set_frequency(float frequency)
  {
    period_ = 1 / frequency;
    omega_ = 2 * PI * frequency;
  }
  
  /**
   * Set the amplitude of the excitation.
   */ 
  inline void set_amplitude(float ampl)
  {
    ampl_ = ampl;
  }

  /**
   * Returns the frequency of the excitation in Hz. 
   */
  inline float get_frequency()
  {
    return 1/period_;
  }

  /**
   * Returns the amplitude of the signal
   */
  inline float get_amplitude()
  {
    return ampl_;
  }

  /**
   * Print a string representation to an ostream.
   */
  virtual ostream& to_string(ostream &os) const;

protected:
  float period_;
  float omega_;
  float ampl_;
};

#endif // EXCITE_EXPSINE_H
