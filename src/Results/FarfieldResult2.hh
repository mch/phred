/* 
   Phred - Phred is a parallel finite difference time domain
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

#ifndef FARFIELD_RESULT2_H
#define FARFIELD_RESULT2_H

#include "Result.hh"
#include "../CSG/CSGBox.hh"

/**
 * This result calculates a near to far field transformation at a
 * number of discreet frequencies and a number of discreet angles. The
 * method implemented here calculates the phasor electric and magnetic
 * currents by taking the DFT of the tangential E and H fields at each
 * time step. At some specified time step, presumably the last, the
 * farfield radiation is calculated. 
 *
 * This does not attempt to normalize to the source. Users will have
 * to do that in post processing for now.
 */ 
class FarfieldResult2 : public Result
{
public:
  FarfieldResult2();
  ~FarfieldResult2();

  /**
   * Allocate memory etc
   */ 
  void init(const Grid &grid);

  /**
   * Deallocate memory etc
   */ 
  void deinit();

  /**
   * Looks at the grid and produces output. 
   *
   * @param grid a reference to a Grid object
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  map<string, Variable *> &get_result(const Grid &grid, 
                                      unsigned int time_step);

  /**
   * Set the start frequency of the range
   */
  inline void set_freq_start(field_t fs)
  {
    freq_start_ = fs;
  }

  /**
   * Set the end frequency of the range
   */
  inline void set_freq_stop(field_t fs)
  {
    freq_stop_ = fs;
  }

  /**
   * Set the number of frequencies of in the range
   */
  inline void set_num_freq(unsigned int nf)
  {
    num_freqs_ = nf;
  }

  /**
   * Get the start frequency of the range
   */
  inline field_t get_freq_start()
  {
    return freq_start_;
  }

  /**
   * Get the end frequency of the range
   */
  inline field_t get_freq_stop()
  {
    return freq_stop_;
  }

  /**
   * Get the number of frequencies of in the range
   */
  inline unsigned int get_num_freq()
  {
    return num_freqs_;
  }

private:
  shared_ptr<CSGBox> box_; /**< The box to use as the surface to
                              integrate currents over. */
  
  region_t grid_box_; /**< The cells that are in our local region */ 

  field_t freq_start_; /**< First frequency in the range */
  field_t freq_stop_; /**< Last frequency in the range */
  unsigned int num_freqs_; /**< Number of frequencies in range */
  field_t freq_space_;

  // Precalculated DFT constants for each frequency

  // Storage space for J and M phasors on each face

  // Storage space for farfield data freq x theta x phi

  // Output variables
  Variable farfield_; 
  Variable freqs_;
  Variable angles_;

};

#endif
