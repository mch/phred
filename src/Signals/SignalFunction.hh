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

#ifndef SIGNAL_H
#define SIGNAL_H

#include "../Types.hh"
#include "../Grid.hh"

/**
 * This class defines a signal_function which can be computed at some
 * timestep on some grid.
 *
 * Subclassess of this are intended to be template parameters to
 * TimeExcitation and SpaceTimeExcitation. 
 */
class Signal
{
private:
protected:
public:
  Signal() 
  {}
  virtual ~Signal() 
  {}

  /**
   * This function is defined in subclasses and produces the signal
   * value to be applied to every point which is excited.
   *
   * @param time the time at which to apply the excitation. Usually
   * just grid.deltat() * time_step.
   */
  virtual field_t signal_function(float time) = 0;

  /**
   * Print a string representation to an ostream.
   */
  virtual ostream& to_string(ostream &os) const
  { return os << "Signal of indeterminate type."; }

  friend ostream& operator<< (ostream& os, const Signal &sf);
};

inline ostream& operator<< (ostream& os, const Signal &sf)
{
  return sf.to_string(os);
}

#endif // SIGNAL_FUNCTION_H
