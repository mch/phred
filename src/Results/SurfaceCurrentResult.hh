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

#ifndef SURFACE_CURRENT_RESULT_H
#define SURFACE_CURRENT_RESULT_H

#include "Result.hh"
#include "../CSG/CSGBox.hh"

/**
 * This class computes the currents on the surface of a box. Time
 * domain and frequency domain results are available. External
 * postprocessors should be able to use this to compute near to far
 * field, or make cool animations.
 */
class SurfaceCurrentResult : public Result
{
public:
  SurfaceCurrentResult();
  ~SurfaceCurrentResult();

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

  field_t *freqs_; /**< Frequency points */

  bool do_dft_; /**< True if we should compute the DFT */ 

  Variable back_Jy_;
  Variable back_Jz_;
  Variable back_My_;
  Variable back_Mz_;

  Variable front_Jy_;
  Variable front_Jz_;
  Variable front_My_;
  Variable front_Mz_;

  Variable left_Jx_;
  Variable left_Jz_;
  Variable left_Mx_;
  Variable left_Mz_;

  Variable right_Jx_;
  Variable right_Jz_;
  Variable right_Mx_;
  Variable right_Mz_;

  Variable bottom_Jx_;
  Variable bottom_Jy_;
  Variable bottom_Mx_;
  Variable bottom_My_;

  Variable top_Jx_;
  Variable top_Jy_;
  Variable top_Mx_;
  Variable top_My_;

  field_t *back_Jy_data_;
  field_t *back_Jz_data_;
  field_t *back_My_data_;
  field_t *back_Mz_data_;

  field_t *front_Jy_data_;
  field_t *front_Jz_data_;
  field_t *front_My_data_;
  field_t *front_Mz_data_;

  field_t *left_Jx_data_;
  field_t *left_Jz_data_;
  field_t *left_Mx_data_;
  field_t *left_Mz_data_;

  field_t *right_Jx_data_;
  field_t *right_Jz_data_;
  field_t *right_Mx_data_;
  field_t *right_Mz_data_;

  field_t *bottom_Jx_data_;
  field_t *bottom_Jy_data_;
  field_t *bottom_Mx_data_;
  field_t *bottom_My_data_;

  field_t *top_Jx_data_;
  field_t *top_Jy_data_;
  field_t *top_Mx_data_;
  field_t *top_My_data_;
  
  
};
