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
class WaveGuideExcitation : public WindowedExcitation
{
public:
  WaveGuideExcitation(SouceFunction *sf);
  ~WaveGuideExcitation();

  virtual field_t window(region_t r, unsigned int x, unsigned int y, 
                         unsigned int z);

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
   * Set the excitation size, location, and orientation. This is
   * totally weak. ugg.
   */
//   void set_wg_excitation(unsigned int x, unsigned int y, unsigned int z,
//                          unsigned int width, unsigned int height, 
//                          Axis a)
//   {
//     x_ = x;
//     y_ = y;
//     z_ = z;
//     width_ = width;
//     height_ = height;
//     axis_ = a;

//     region_.xmin = x_;
//     region_.ymin = y_;
//     region_.zmin = z_;

//     switch (axis_)
//     {
//     case X_AXIS:
//       region_.ymax = y_ + width;
//       region_.zmax = z_ + height;
//       break;
//     case Y_AXIS:
//       region_.xmax = x_ + width;
//       region_.zmax = z_ + height;
//       break;
//     case Z_AXIS:
//       region_.xmax = x_ + width;
//       region_.ymax = y_ + height;
//       break;
//     }
//   }



  /**
   * From the window size, placement, and orientation, construct a
   * suitable region and set the mode numbers for each axis. Make
   * sure everything is inside the local region, etc.
   */ 
  //virtual void init(const Grid &grid);

protected:
  unsigned int mode_x_;
  unsigned int mode_y_;
  unsigned int mode_z_;

  //unsigned int x_, y_, z_;
  //unsigned int width_, height_;
  //Axis axis_;

private:
};

#endif WAVEGUIDE_EXCITATION_H
