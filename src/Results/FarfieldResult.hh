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
 * This class implements Luebbers' method, Taflove, Computational
 * Electrodynamics, 2nd ed, pg 366
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

  // Number of time steps it takes for a wave to progagate through
  // freespace between the two points which are farthest from each
  // other on the Huygen's surface.
  unsigned int t_cross_;

  // Field components of W and U. These are each 3 dimensional grids,
  // indexed by phi, theta, and farfield tstep. These are contiguous
  // along the time step dimension.
  field_t *Wx_, *Wy_, *Wz_;
  field_t *Ux_, *Uy_, *Uz_;

  // When updating W and U, it is necessary to approximate the
  // derivative of E and H w.r.t. time. For this, it is necessary to
  // store the past values of E and H at each point on the surface. (I
  // think... maybe theres something cleaver that can be done to avoid
  // this....)
  field_t *E_temp_, *H_temp_;

  // E_theta and E_phi values, calculated by get_post_result(). These
  // are 3d, indexed by phi, theta, and finally contigouos along
  // fftimestep.
  field_t *E_theta_, *E_phi_;

  // Output variables
  Variable freqs_;
  Variable theta_;
  Variable phi_;
  Variable E_theta_var_;
  Variable E_phi_var_;

  /**
   * A helper for accessing data in the E and H temporary arrays. Data
   * is treated as contiguous along x2.
   * 
   * @param face the face, 0-5
   * @param component the field component, x=0, y=1, z=2
   * @param x1 first coordinate on the face, beginning at zero
   * @param x2 second coordinate on the face, beginning at zero
   * @param size1 length of the first dimension of the face
   * @param size2 length of the second dimension of the face
   */
  static inline int temp_index(int face, int component, 
                               int x1, int x2,
                               int size1, int size2)
  {
    return x2 + (x1 + (component + face * 3) * size1) * size2;
  }

  /**
   * A helper for accessing data in the W and U arrays. Data is
   * contiguous along the time step.
   *
   * @param phi observation point angle with the z axis
   * @param theta observation point angle with the x axis
   * @param fftstep far field time step
   */ 
  inline int WU_index(int phi_idx, int theta_idx, int ffstep)
  {
    return ffstep + (theta_idx + phi_idx * theta_data_.length()) * ff_tsteps_;
  }

};

#endif // FARFIELD_RESULT_H
