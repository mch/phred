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
 *
 * \bug DFT is not yet implemented. 
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
   * Set the CSGBox to calculate the currents on
   */
  inline void set_region(shared_ptr<CSGBox> box)
  {
    box_ = box;
  }

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

  /**
   * Print a string representation to an ostream.
   */
  virtual ostream& to_string(ostream &os) const;

private:
  shared_ptr<CSGBox> box_; /**< The box to use as the surface to
                              integrate currents over. */
  
  shared_ptr<Block> region_; /**< The cells that are in our local region */ 

  field_t freq_start_; /**< First frequency in the range */
  field_t freq_stop_; /**< Last frequency in the range */
  unsigned int num_freqs_; /**< Number of frequencies in range */
  field_t freq_space_;

  field_t *freqs_; /**< Frequency points */

  bool do_dft_; /**< True if we should compute the DFT */ 

  Variable Jt1_[6];
  Variable Jt2_[6];
  Variable Mt1_[6];
  Variable Mt2_[6];

  field_t *Jt1_data_[6];
  field_t *Jt2_data_[6];
  field_t *Mt1_data_[6];
  field_t *Mt2_data_[6];


  // Templated function that can use adapter classes
  template<class T>
  inline void calc_currents(unsigned int xmin, unsigned int xmax, 
                            unsigned int ymin, unsigned int ymax, 
                            unsigned int zmin, unsigned int zmax, 
                            int face_idx, const Grid &grid);  
};

#endif
