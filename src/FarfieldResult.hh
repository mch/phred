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

#ifndef FARFIELD_RESULT_H
#define FARFIELD_RESULT_H

#include "Result.hh"

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
 * Almost all of this code (except for the MPI stuff) is borrowed from
 * Jan's implementation.
 *
 * \bug Support for more than one rank (parallel communication) is
 * not yet implemented
 *
 * \bug Current implementation is pretty hard on memory... mostly
 * because of the need to copy things into a contigous output array.
 */
class FarfieldResult : public Result
{
private:
protected:
  float theta_start_;
  float theta_stop_;
  float phi_start_;
  float phi_stop_;
  unsigned int num_pts_; /**< Number of points on arc */

  Axis axis_; /**< Axis to rotate around if it can't be deduced from
                 theta and phi start/stop */

  field_t freq_start_; /**< First frequency in the range */
  field_t freq_stop_; /**< Last frequency in the range */
  unsigned int num_freqs_; /**< Number of frequencies in range */
  field_t freq_space_;
  
  region_t global_r_; /**< Region specifying surface of Huygen's box in
                         the global grid (this result only deals in
                         the global grid) */

  float *theta_;
  float *phi_;

  float **e_theta_re_;
  float **e_theta_im_;
  float **e_phi_re_;
  float **e_phi_im_;
  float ***jff_mom_;
  float ***mff_mom_;

  int t_cross_; /**< ?? */
  
  int rank_; /**< The rank of this process */
  int size_; /**< Number of ranks in the MPI communicator */

  FfType output_type_; /**< Type of result to produce */

  field_t result_; /**< Data to output */

  /**
   * Helper function to calculate the result. 
   */
  void ffpu_calculate(const Grid &grid, unsigned int time_step);

  /**
   * Calculate the updated moments. 
   */
  void ffpu_moment_update(const Grid &grid);

  /**
   * I don't know what this one does!
   */
  void ffpu_ntff(float theta, float phi, float *ms, float *js,
                 float *e_theta, float *e_phi);

  /**
   * I don't know what this one does!
   */
  void ffpu_moment_cycle();

public:
  FarfieldResult();
  ~FarfieldResult();

  /**
   * Set the type of result to produce. Defaults to ETHETAPHI.
   */
  inline void set_output_type(FfType ot)
  {
    output_type_ = ot;
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
   * Set the rank and size of the MPI communicator. This is required
   * for this Resault, because it needs to fetch data from other
   * ranks to do the computation. 
   */ 
  void set_mpi_rank_size(int rank, int size);

  /**
   * Set the angles that define the arc on which the transformation
   * will be calculated. If the arc cannot be uniquely determined
   * from these angles, also set the axis of rotation (which defaults
   * to X_AXIS). 
   */ 
  inline void set_angles(float theta_start, float theta_stop, 
                         float phi_start, float phi_stop, 
                         unsigned int num_pts)
  {
    theta_start_ = theta_start;
    theta_stop_ = theta_stop;
    phi_start_ = phi_start;
    phi_stop_ = phi_stop;
    num_pts_ = num_pts;
  }

  /**
   * Set the axis of rotation, which is only used if the arc cannot
   * be uniquely determined from the angles. Defaults to X_AXIS. 
   */
  inline void set_axis(Axis a)
  {
    axis_ = a;
  }

  /**
   * Compute the near to farfield transformation. 
   *
   * @param grid a reference to a Grid object
   * @param time_step current time step
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  Data &get_result(const Grid &grid, unsigned int time_step);

  /**
   * Setup the result, allocate memory, etc. Called just before the
   * simulation starts. 
   */
  void init(const Grid &grid);
  
  /**
   * Deallocates memory. 
   */
  void deinit(const Grid &grid);

protected:
  /**
   * Computes arc angles
   */ 
  void arc_connect();

  /**
   * Determines if the theta and phi values lead to a well defined
   * arc. If not, the axis of rotation must be relied upon.
   */
  int is_parallel(float th1, float ph1, float th2, float ph2);

  /**
   * Avoids singularities in atan.
   */
  float atan2_local(float y, float x);

};

#endif // FARFIELD_RESULT_H
