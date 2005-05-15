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

#ifndef PERIODIC_EXCITE_H
#define PERIODIC_EXCITE_H

#include "Excitation.hh"

/**
 * This class defines an excitation to be used by periodic boundary
 * conditions. 
 */ 
class PeriodicExcitation : public Excitation {
public:

  PeriodicExcitation(shared_ptr<Signal> sf);

  virtual ~PeriodicExcitation();

  /**
   * This function is called to apply this excitation to the grid. 
   *
   * @param grid the grid to operate on
   * @param time_step the time step we are on
   * @param type the field type we are allowed to excite. 
   */
  void excite(Grid &grid, unsigned int time_step, 
              FieldType type);

  /**
   * Ensure that the excitation is being applied to an entire plane. 
   */
  void init(const Grid &grid);

  /**
   * Print a string representation to an ostream.
   */
  ostream& to_string(ostream &os) const;

  /**
   * Set the CSGBox inside which the excitation is to be applied. 
   *
   * @param box the CSGBox to apply the excitation to. 
   */
  void set_region(shared_ptr<CSGBox> box, Face face);

  /**
   * Returns the face
   */ 
  Face get_face() const
  { return face_; }

private:
  void set_region(shared_ptr<CSGBox> box) {}

  Face face_; /**< The face to apply the excitation to */ 
};

#endif // PERIODIC_EXCITE_H
