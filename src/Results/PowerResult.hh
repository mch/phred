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

#ifndef POWER_RESULT_H
#define POWER_RESULT_H

#include <complex>

#include "DFTResult.hh"
#include "../GridPlane.hh"

/**
 * Calculates the power flowing through a surface at specific
 * frequencies. 
 *
 * \bug This now uses the averaging functions of the GridPlane. The
 * region the power is taken over must allow room for the averaging
 * functions to work. 
 */
class PowerResult : public DFTResult
{
public:

  PowerResult();

  /**
   * Constructor. Sets the frequency parameters at the same time. 
   *
   * @param freq_start first frequency in the range
   * @param freq_stop last frequency in the range
   * @param num_freqs total number of frequencies to report data for
   */
  PowerResult(field_t freq_start, field_t freq_stop, 
              unsigned int num_freqs);
  ~PowerResult();

  /**
   * Allocate memory, sanity checks. 
   */ 
  virtual void init(const Grid &grid);
  
  /**
   * Free memory etc
   */
  virtual void deinit();

  /**
   * Looks at the grid and produces output
   *
   * @param grid a reference to a Grid object
   */
  void calculate_result(const Grid &grid, 
                        unsigned int time_step);

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
  field_t time_power_; /**< Power at the current instant in time domain */ 

  shared_ptr<CSGBox> box_; 
  Face face_;
  shared_ptr<Block> region_;
  bool has_data_; // True if this local region generates data. 

  field_t cell_area_; /**< Area for each FDTD cell... */ 

  Variable real_var_; /**< Power data */
  Variable imag_var_; /**< Power data */
  Variable freq_var_; /**< Frequency data */
  Variable power_var_; /**< Power at a time instant */ 

//   field_t *et1r_;
//   field_t *et1i_;
//   field_t *ht1r_;
//   field_t *ht1i_;

//   field_t *et2r_;
//   field_t *et2i_;
//   field_t *ht2r_;
//   field_t *ht2i_;

  complex<field_t> *et1_;
  complex<field_t> *et2_;
  complex<field_t> *ht1_;
  complex<field_t> *ht2_;

  unsigned int x_size_;
  unsigned int y_size_;
  unsigned int z_size_;

  // Need to time average E field to get the right power.
  field_t *prev_et1_; /**< Previous value of 1st tangential E component */ 
  field_t *prev_et2_; /**< Previous value of 2nd tangential E component */ 

  /**
   * Print a string representation to an ostream.
   */
  virtual ostream& to_string(ostream &os) const;

private:
};

#endif
