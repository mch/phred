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

/** \class Gaussm
 * \brief A modulated gaussian excitation
 *
 * A Gaussian function modulated by a sine wave. 
 */

#ifndef EXCITE_GAUSS_PULSE_H
#define EXCITE_GAUSS_PULSE_H

#include "Signal.hh"
#include "../Constants.hh"

#include <cmath>

class GaussPulse : public Signal
{
protected:
  field_t alpha_;
  field_t tau_p_;
  
public:
  GaussPulse();
  virtual ~GaussPulse();

  /**
   * Set the parameters of the Gauss pulse.
   *
   * @param alpha scale factor (defaults to 1)
   * @param tau_p characteristic time
   */
  void set_parameters(field_t alpha, field_t tau_p);

  /**
   * Get the current alpha value
   *
   * @return alpha
   */
  field_t get_alpha() const;

  /**
   * Get the current characteristic time
   *
   * @return characteristic time
   */
  field_t get_taup() const;

  /**
   * Produces a modulated gauss function. 
   *
   * @param time the time at which to apply the excitation
   *
   * @return the value of the excitation
   */
  field_t signal_function(float time) const;

  /**
   * Returns the amount of time it will take for the Gaussian to rise
   * and fall.
   */ 
  field_t length() const;

  /**
   * Print a string representation to an ostream.
   */
  virtual ostream& to_string(ostream &os) const;

};

#endif // EXCITE_GAUSSM_H
