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

#ifndef WAVEGUIDE_EXCITATION_H
#define WAVEGUIDE_EXCITATION_H

#include "WindowedExcitation.hh"

/**
 * A wave guide excitation. Exciates TE_{nm} modes across a retangular
 * window.
 *
 * \bug this class really isn't that intellegent. In fact, it's
 * downright stupid.
 */ 
class WaveguideExcitation : public WindowedExcitation
{
public:
  WaveguideExcitation(shared_ptr<SourceFunction> sf);
  ~WaveguideExcitation();

  virtual field_t window(float x, float y, float z);

  /**
   * Set the mode to excite. 
   *
   * @param n the number of half wavelengths in the x direction
   * @param m the number of half wavelengths in the y direction
   * @param o the number of half wavelengths in the z direction
   */
  void set_mode(unsigned int n, unsigned int m, unsigned int o)
  {
    mode_x_ = n;
    mode_y_ = m;
    mode_z_ = o;
  }

  /**
   * Print a string representation to an ostream.
   */
  virtual ostream& to_string(ostream &os) const;

protected:
  unsigned int mode_x_;
  unsigned int mode_y_;
  unsigned int mode_z_;

private:
};

#endif // WAVEGUIDE_EXCITATION_H
