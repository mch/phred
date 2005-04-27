/* 
   Phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004-2005 Matt Hughes <mhughe@uvic.ca>

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

#ifndef CHECKPOINT_H
#define CHECKPOINT_H

/**
 * This class is an ABC which must be implemented by any class wishing
 * to checkpoint its state. This is required because some computer
 * systems do not allow programs to run for longer than a preset
 * limit, such as 24 hours. This also allows a job to be stopped and
 * restarted at a later time. It can also be used to recover from a
 * crash.
 *
 * This class is intedent to be subclasses as part of multiple
 * inheritance, in particular with the LifeCycle class. 
 *
 * \bug THIS IS JUST A SHELL for now... flesh this out.
 */ 
class Checkpoint
{
public:
  virtual ~Checkpoint()
  {}

  /**
   * Subclasses must implement this to save data...
   *
   * @param ...
   */ 
  void checkpoint_save()
  {}

  /**
   * Subclasses must implement this to load data and resume from a
   * checkpoint.
   *
   * @param ...
   */ 
  void checkpoint_load()
  {}

private:

};
