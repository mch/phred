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

#ifndef SIGNAL_TIME_RESULT_H
#define SIGNAL_TIME_RESULT_H

#include "Result.hh"
#include "../Signals/Signal.hh"

/**
 * Outputs a signal function at the timestep.
 */
class SignalTimeResult : public Result
{
private:
protected:
  Signal &te_; /**< The time excitation to save. */

  field_t result_[2]; /**< Storage for the result. */
  
  Variable var_;

public:
  SignalTimeResult(Signal &te);
  ~SignalTimeResult();

  /**
   * Set the name on the variable
   */
  virtual void init(const Grid &grid);

  /**
   * Set the signal or Signal to use.
   * @param a reference to a Signal object
   */ 
  void set_excitation(const Signal &te);

  /**
   * Produces output from the signal, considering only the time step. 
   *
   * @param grid a reference to a Grid object
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  virtual map<string, Variable *> &get_result(const Grid &grid, 
                                              unsigned int time_step);

  /**
   * Print a string representation to an ostream.
   */
  virtual ostream& to_string(ostream &os) const;

};

#endif // SIGNAL_TIME_RESULT_H
