/* 
   phred - Phred is a parallel finite difference time domain
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

#ifndef FAKE_RESULT_H
#define FAKE_RESULT_H

#include "Result.hh"

/**
 * Just for testing
 */ 
class FakeResult : public Result
{
public:
  FakeResult();
  ~FakeResult();

  void init(const Grid &grid);
  void deinit();

  /**
   * Print a string representation to an ostream.
   */
  virtual ostream& to_string(ostream &os) const;

  void calculate_result(const Grid &grid, unsigned int time_step);

protected:
  Variable var_;
  field_t *data_;
  MPI_Datatype dtype_;

};

#endif
