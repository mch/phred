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

#include "../config.h"
#include "DFTResult.hh"
#include "../CSG/CSGBox.hh"

#ifdef HAVE_COMPLEX
#include <complex>
#endif

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
class FarfieldResult2 : public DFTResult
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
   * Set the range of theta (angle from the X axis) to calculate the
   * farfield data for. The range must not span more than 360 degrees. 
   *
   * @param theta_start starting angle in radians
   * @param theta_stop end angle in radians
   * @param num_theta number of points of theta to calculate FF data
   * for. Must be greater than 0.
   */ 
  void set_theta(field_t theta_start, field_t theta_stop, 
                 unsigned int num_theta);

  /**
   * Set the range of phi (angle from the Z axis) to calculate the
   * farfield data for. The range must not span more than 180 degrees. 
   *
   * @param phi_start starting angle in radians
   * @param phi_stop end angle in radians
   * @param num_phi number of points of phi to calculate FF data
   * for. Must be greater than 0.
   */ 
  void set_phi(field_t phi_start, field_t phi_stop, 
               unsigned int num_phi);

  /**
   * Set the range of theta (angle from the X axis) to calculate the
   * farfield data for. The range must not span more than 360 degrees. 
   *
   * @param theta_start starting angle in radians
   * @param theta_stop end angle in radians
   * @param num_theta number of points of theta to calculate FF data
   * for. Must be greater than 0.
   */ 
  void set_theta_degrees(field_t theta_start, field_t theta_stop, 
                         unsigned int num_theta);

  /**
   * Set the range of phi (angle from the Z axis) to calculate the
   * farfield data for. The range must not span more than 180 degrees. 
   *
   * @param phi_start starting angle in degrees
   * @param phi_stop end angle in degrees
   * @param num_phi number of points of phi to calculate FF data
   * for. Must be greater than 0.
   */ 
  void set_phi_degrees(field_t phi_start, field_t phi_stop, 
                       unsigned int num_phi);

  /**
   * Print a string representation to an ostream.
   */
  ostream& to_string(ostream &os) const;

private:
  shared_ptr<CSGBox> box_; /**< The box to use as the surface to
                              integrate currents over. */
  
  region_t grid_box_; /**< The cells that are in our local region */ 

  // Angle range
  Interval<field_t> theta_data_;
  Interval<field_t> phi_data_;

  // Storage space for J and M phasors on each face
#ifdef HAVE_COMPLEX
  complex<field_t> *Jt1_data_[6];
  complex<field_t> *Jt2_data_[6];

  complex<field_t> *Mt1_data_[6];
  complex<field_t> *Mt2_data_[6];

  // Storage space for farfield E_theta data theta x phi x freq
  complex<field_t> *e_theta_data_;

  // Storage space for farfield E_phi data theta x phi x freq
  complex<field_t> *e_phi_data_;

  // Storage space for farfield H_theta data theta x phi x freq
  complex<field_t> *h_theta_data_;

  // Storage space for farfield H_phi data theta x phi x freq
  complex<field_t> *h_phi_data_;
#endif

  // Storage space for farfield RCS data theta x phi x freq
  field_t *rcs_data_;

  // Output variables
  Variable rcs_;
//   Variable E_phi_;
//   Variable E_theta_;
//   Variable H_phi_;
//   Variable H_theta_;
  Variable freqs_;
  Variable theta_;
  Variable phi_;

};

#endif
