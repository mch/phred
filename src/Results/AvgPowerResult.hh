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

#ifndef AVG_POWER_RESULT_H
#define AVG_POWER_RESULT_H

#include <complex>

#include "DFTResult.hh"
#include "../GridPlane.hh"

/**
 * Calculates the average power flowing through a surface at specific
 * frequencies. 
 */
class AvgPowerResult : public DFTResult
{
public:

  AvgPowerResult();

  /**
   * Constructor. Sets the frequency parameters at the same time. 
   *
   * @param freq_start first frequency in the range
   * @param freq_stop last frequency in the range
   * @param num_freqs total number of frequencies to report data for
   */
  AvgPowerResult(field_t freq_start, field_t freq_stop, 
              unsigned int num_freqs);
  ~AvgPowerResult();

  /**
   * Allocate memory, sanity checks. 
   */ 
  void init(const Grid &grid);
  
  /**
   * Free memory etc
   */
  void deinit();

  /**
   * Adds DFT'd values for this time step
   *
   * @param grid a reference to a Grid object
   */
  void calculate_result(const Grid &grid, 
                        unsigned int time_step);

  /**
   * Calculates the output average S from the stored DFT'ed E and H fields. 
   */
  void calculate_post_result(const Grid &grid);

  /**
   * Set the plane to do the calculation on
   */
  inline void set_region(shared_ptr<CSGBox> box, Face face)
  {
    face_ = face;
    box_ = box;
  }

protected:
  field_t *power_real_; /**< Power at each frequency */ 
  field_t *power_imag_; /**< Power at each frequency */ 

  complex<field_t> *et1_;
  complex<field_t> *et2_;
  complex<field_t> *ht1_;
  complex<field_t> *ht2_;

  shared_ptr<CSGBox> box_; 
  Face face_;
  shared_ptr<Block> region_;
  bool has_data_; // True if this local region generates data. 

  field_t cell_area_; /**< Area for each FDTD cell... */ 

  Variable real_var_; /**< Power data */
  Variable imag_var_; /**< Power data */
  Variable freq_var_; /**< Frequency data */

  unsigned int x_size_;
  unsigned int y_size_;
  unsigned int z_size_;

  // Need to time average E field to get the right power.
  field_t *prev_et1_; /**< Previous value of 1st tangential E component */ 
  field_t *prev_et2_; /**< Previous value of 2nd tangential E component */ 

  // Storage for pre-calculated sin and cos values
  field_t *cos_temp_;
  field_t *sin_temp_;
  
  field_t *e_cos_temp_;
  field_t *e_sin_temp_;

  /**
   * Print a string representation to an ostream.
   */
  ostream& to_string(ostream &os) const;

private:
};

#endif
