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

#ifndef SOURCE_TIME_RESULT_H
#define SOURCE_TIME_RESULT_H

#include "Result.hh"
#include "SourceFunction.hh"

/**
 * Outputs a source function at the timestep.
 */
class SourceTimeResult : public Result
{
private:
protected:
  SourceFunction &te_; /**< The time excitation to save. */

  field_t result_[2]; /**< Storage for the result. */
  
  Variable var_;

public:
  SourceTimeResult(SourceFunction &te);
  ~SourceTimeResult();

  /**
   * Set the source or SourceFunction to use.
   * @param a reference to a SourceFunction object
   */ 
  void set_excitation(const SourceFunction &te);

  /**
   * Produces output from the source, considering only the time step. 
   *
   * @param grid a reference to a Grid object
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  virtual map<string, Variable *> &get_result(const Grid &grid, 
                                              unsigned int time_step);

};

#endif // SOURCE_TIME_RESULT_H
