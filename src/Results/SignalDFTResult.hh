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

#ifndef SIGNAL_DFT_RESULT_H
#define SIGNAL_DFT_RESULT_H

#include "DFTResult.hh"
#include "../Signals/Signal.hh"

/**
 * Outputs the DFT of a signal function. Only really applies to
 * excitations that are applied at a single point in space.
 */
class SignalDFTResult : public DFTResult
{
private:
protected:
  Signal &te_; /**< The dft excitation to save. */

  field_t *result_; /**< Storage for the result. Interlevaved data;
                       freq, Real DFT value, Imag DFT value, etc */

  Variable var_;

public:
  SignalDFTResult(Signal &te);
  SignalDFTResult(Signal &te, field_t freq_start,
                  field_t freq_stop, unsigned int num_freqs);
  ~SignalDFTResult();

  /**
   * Set the signal or TimeExcitation to use.
   * @param a reference to a TimeExcitation object
   */ 
  void set_excitation(const Signal &te);

  /**
   * Produces a running DFT from the signal, considering only the time
   * step. Output is only available at time_end. 
   *
   * @param grid a reference to a Grid object
   * @param time_step current time step
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  virtual map<string, Variable *> &get_result(const Grid &grid, 
                                              unsigned int time_step);

  /**
   * Setup the result, allocate memory, etc. Called just before the
   * simulation starts. 
   */
  void init(const Grid &grid);
  
  /**
   * Deallocates memory. 
   */
  void deinit();

  /**
   * Print a string representation to an ostream.
   */
  virtual ostream& to_string(ostream &os) const;

};

#endif // SIGNAL_DFT_RESULT_H
