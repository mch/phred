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

#ifndef WINDOWED_EXCITATION_H
#define WINDOWED_EXCITATION_H

#include "Excitation.hh"
#include "../Globals.hh"

/**
 * An excitation that is windowed in 3 space by some windowing function. 
 */
class WindowedExcitation : public Excitation
{
private:
protected:

  // Box ranges in real global coordinates. For help in calculating
  // the window.
  float xmin_, xmax_, ymin_, ymax_, zmin_, zmax_;

  // The part of the box that is in the local domain. For looping
  // across the cells contained in the box.
  float lxmin_, lymin_, lzmin_;

public:
  WindowedExcitation(shared_ptr<SignalFunction> sf);

  virtual ~WindowedExcitation();

  /**
   * init the windowing function.
   */ 
  virtual void init(const Grid &grid);

  /**
   * Subclasses must override this to provide the window
   * function. This is a bit on the wasty side, because it
   * unnecessarilly calculates the x and y parts of the window at
   * every z step. If profiling shows there is too much time wasted
   * here, then optimize this by spliting it into three functions. 
   *
   * @param x the x position in the global real coordinate system
   * @param y the y position in the global real coordinate system
   * @param z the z position in the global real coordinate system
   * @return a value describing the strength of the field at the given point. 
   */
  virtual field_t window(float x, float y, float z) = 0;

  /**
   * This applied the windowed excitation to the grid.
   *
   * @param grid the grid to mess with
   * @param time_step the time step we are on
   * @param type the field components to excite (E or H)
   */
  virtual void excite(Grid &grid, unsigned int time_step, 
                      FieldType type);
  
};

#endif // WINDOWED_EXCITATION_H
