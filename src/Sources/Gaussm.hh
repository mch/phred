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

/** \class Gaussm
 * \brief A modulated gaussian excitation
 *
 * A Gaussian function modulated by a sine wave. 
 */

#ifndef EXCITE_GAUSSM_H
#define EXCITE_GAUSSM_H

#include "SourceFunction.hh"
#include "../Constants.hh"

#include <math.h>

class Gaussm : public SourceFunction
{
protected:
  field_t alpha_;
  field_t deltaf_;
  field_t f0_;
  
public:
  Gaussm();
  virtual ~Gaussm();

  /**
   * Set the parameters of the modulated Gauss function.
   *
   * @param alpha scale factor (defaults to 1)
   * @param deltaf the frequency width of the pulse (defaults to 50 MHz)
   * @param f0 the centre frequency of the pulse (defaults to 1 GHz)
   */
  void set_parameters(field_t alpha, field_t deltaf, field_t f0);

  /**
   * Get the current alpha value
   *
   * @return alpha
   */
  field_t get_alpha();

  /**
   * Get the current frequency range
   *
   * @return frequency range
   */
  field_t get_deltaf();

  /**
   * Get the current centre frequency
   *
   * @return centre frequency
   */
  field_t get_f0();

  /**
   * Produces a modulated gauss function. 
   *
   * @param time the time at which to apply the excitation
   *
   * @return the value of the excitation
   */
  field_t source_function(float time);

};

#endif // EXCITE_GAUSSM_H
