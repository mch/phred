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

#ifndef LIFE_CYCLE_H
#define LIFE_CYCLE_H

class Grid;

/**
 * An interface that defines life cycle functions used by a FDTD sim
 * management object to set up member objects when leaving define mode
 * so that they can allocate memory and init constants that may depend
 * on the state of the grid or other objects. 
 */
class LifeCycle
{
private:
protected:
public:
  LifeCycle();

  virtual ~LifeCycle() = 0;

  /**
   * Subclasses can implement this to allocate memory or init
   * constants or whatever. Called just before the simulation
   * starts. The default implementation does nothing. 
   */
  virtual void init(const Grid &grid)
  {}
  
  /**
   * Subclasses can implement this to deallocate memory or
   * whatever. Called just after the simulation ends. The default
   * implementation does nothing. 
   */
  virtual void deinit()
  {}
};

#endif // LIFE_CYCLE_H
