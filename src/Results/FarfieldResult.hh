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

#ifndef FARFIELD_RESULT_H
#define FARFIELD_RESULT_H

#include "DFTResult.hh"
#include "../config.h"

/**
 * Available output types
 */
enum FfType {
  ETHETA,
  EPHI,
  HTHETA, 
  HPHI, 
  ETHETAPHI, 
  HTHETAPHI, 
  RCSNORM, 
  RCSDBPOL
};

/**
 * Computes near field to far field transformation. This result must
 * collect data from other ranks using MPI communication. This might
 * be refactored later...
 *
 * This class implements Luebbers' method, R. J. Luebbers,
 * M. Scheider, "A finite-difference time-domain near zone to far zone
 * transformation", IEEE Transactions on Antennas and Propagation, Vol
 * 39, No 4, April 1991, 429-433.
 */
class FarfieldResult : public DFTResult
{
public:
  FarfieldResult();
  ~FarfieldResult();

  /**
   * Compute the near to farfield transformation. 
   *
   * @param grid a reference to a Grid object
   * @param time_step current time step
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  map<string, Variable *> &get_result(const Grid &grid, 
                                      unsigned int time_step);

  /**
   * Returns the frequency range, theta, and phi ranges. 
   */ 
  map<string, Variable *> &get_pre_result(const Grid &grid);

  /**
   * Returns the farfield RCS data. 
   */ 
  map<string, Variable *> &get_post_result(const Grid &grid);

  /**
   * Set the box over which the currents should be calculated so that
   * the farfield can be calculated from those currents.
   */ 
  void set_region(shared_ptr<CSGBox> box)
  { box_ = box; }

  /**
   * Set the range of theta (angle from the Z axis) to calculate the
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
   * Set the range of phi (angle from the X axis) to calculate the
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
   * Set the range of theta (angle from the Z axis) to calculate the
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
   * Set the range of phi (angle from the X axis) to calculate the
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
   * Setup the result, allocate memory, etc. Called just before the
   * simulation starts. 
   */
  void init(const Grid &grid);
  
  /**
   * Deallocates memory. 
   */
  void deinit();

  /**
   * Set the radius of the observation sphere
   */ 
  void set_radius(field_t r);

  /**
   * Returns the radius of the observation sphere
   */ 
  inline field_t get_radius()
  { return r_; }

  /**
   * Print a string representation to an ostream.
   */
  virtual ostream& to_string(ostream &os) const;

protected:
  shared_ptr<CSGBox> box_; /**< The box to use as the surface to
                              integrate currents over. */
  
  shared_ptr<Block> region_; /**< Grid cells in the local grid */ 

  // Angle range
  Interval<field_t> theta_data_;
  Interval<field_t> phi_data_;

  // Distance to observation point... hmm...
  field_t r_;

  // Number of farfield timesteps we need to save data for
  unsigned int ff_tsteps_;

  // Time it takes for a wave to progagate through freespace between
  // the two points which are farthest from each other on the Huygen's
  // surface.
  field_t t_cross_;

  // Field components of W and U. These are each 3 dimensional grids,
  // indexed by phi, theta, and farfield tstep. These are contiguous
  // along the time step dimension.
  field_t *Wx_, *Wy_, *Wz_;
  field_t *Ux_, *Uy_, *Uz_;

  // E_theta and E_phi values, calculated by get_post_result(). These
  // are 3d, indexed by phi, theta, and finally contigouos along
  // fftimestep.
  field_t *E_theta_, *E_phi_;

  // DFT'ed normalized farfield E_theta and E_phi. These are 3d,
  // indexed by phi, theta, and finally contigouos along frequency. 
  field_t *e_theta_real_, *e_theta_imag_;
  field_t *e_phi_real_, *e_phi_imag_;

  // DFT'ed farfield RCS, 3d, as above. 
  field_t *rcs_;

  // Output variables
  Variable freqs_;
  Variable theta_;
  Variable phi_;
  Variable E_theta_var_;
  Variable E_phi_var_;

  // DFT'd variables
  Variable e_tr_;
  Variable e_ti_;
  Variable e_pr_;
  Variable e_pi_;
  Variable rcs_var_;
  
  void idx_tests();
  void dump_temps();

  /**
   * A helper for calculating the offset into the dft arrays
   *
   * @param phi observation point angle with the x axis
   * @param theta observation point angle with the z axis
   * @param f_idx frequency index
   */
  inline int dft_index(int phi_idx, int theta_idx, int f_idx)
  {
    return f_idx + (theta_idx + phi_idx * theta_data_.length()) 
      * frequencies_.length();
  }

};

/**
 * A helper for accessing data in the W and U arrays. Data is
 * contiguous along the time step.
 *
 * @param phi observation point angle with the x axis
 * @param theta observation point angle with the z axis
 * @param fftstep far field time step
 */ 
inline static int WU_index(int phi_idx, int theta_idx, int ffstep, 
                           int theta_length, int ff_tsteps)
{
  return ffstep + (theta_idx + phi_idx * theta_length) * ff_tsteps;
}


#endif // FARFIELD_RESULT_H
