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

#ifndef TOTAL_SCATTERED_H
#define TOTAL_SCATTERED_H

#include "Excitation.hh"

/**
 * This class implements the total/scattered field formulation as
 * described by Taflove. The source function is uniformly applied at
 * one face. Outgoing fields have the incident field subtracted from
 * them, so that only the scattered fields leave this box.
 */ 
class TotalScattered : public Excitation
{
public:
  TotalScattered(shared_ptr<Signal> sf);
  ~TotalScattered();

  /**
   * This function applies the total / scattered field formulation to
   * the grid. It excites on one face, updates the stored incident
   * field time steps, and subtracts the incident field from the
   * scattered field at the outer boundaries.
   *
   * @param grid the grid to operate on
   * @param time_step the time step we are on
   * @param type the field type we are allowed to excite. 
   */
  void excite(Grid &grid, unsigned int time_step, 
              FieldType type);

  /**
   * Initialize the memory required to timestep the initial field
   * outside of the grid, etc.
   */ 
  void init(const Grid &grid);

  /**
   * Deinit memory. 
   */ 
  void deinit();

  /**
   * Print a string representation to an ostream.
   */
  virtual ostream& to_string(ostream &os) const;

private:

};

#endif // TOTAL_SCATTERED_H
