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

#ifndef SOURCE_FUNCTION_H
#define SOURCE_FUNCTION_H

#include "../Types.hh"
#include "../Grid.hh"

/**
 * This class defines a source_function which can be computed at some
 * timestep on some grid.
 *
 * Subclassess of this are intended to be template parameters to
 * TimeExcitation and SpaceTimeExcitation. 
 */
class SourceFunction
{
private:
protected:
public:
  SourceFunction() 
  {}
  virtual ~SourceFunction() 
  {}

  /**
   * This function is defined in subclasses and produces the source
   * value to be applied to every point which is excited.
   *
   * @param grid the grid to which the function is being applied. Only for
   * things such as delta_t, the result of this function is applied in
   * the excite() function (in order to factor out the loop). 
   *
   * @param time_step the time at which to apply the excitation
   */
  virtual field_t source_function(const Grid &grid, float time_step) = 0;

};

#endif // SOURCE_FUNCTION_H
