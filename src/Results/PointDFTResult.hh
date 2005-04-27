/* 
   phred - Phred is a parallel finite difference time domain
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

#ifndef POINT_DFT_RESULT_H
#define POINT_DFT_RESULT_H

#include "DFTResult.hh"

/**
 * Outputs the DFT of the a field component at a point in space. 
 */
class PointDFTResult : public DFTResult
{
private:
protected:
  field_t *result_; /**< Storage for the result. Interlevaved data;
                       freq, Real DFT value, Imag DFT value, etc */

  field_t prev_e_[3]; /**< Previous values of the E components,
                         because it is necessary to average over two
                         time steps to get the value at dt*tstep,
                         where the H components are.  */

  point space_point_; /**< The point in real space to aquire data at. */ 

  grid_point point_; /**< The point to aquire the data at (global) */
  grid_point l_; /**< The point to aquire the data at (local) */

  bool ours_; /**< True if the point is in our part of the grid. */

  Variable var_; /**< Our variable */

public:
  PointDFTResult();
  PointDFTResult(field_t freq_start, field_t freq_stop, 
                 unsigned int num_freqs);
  ~PointDFTResult();

  /**
   * Set the point in global coordinates
   *
   * @param point
   */
  inline void set_point(point p)
  {
    space_point_ = p;
  }

  /**
   * Get the point in global coordinates
   *
   * @return grid_point
   */
  inline point get_point()
  {
    return space_point_;
  }

  /**
   * Produces a running DFT from the source, considering only the time
   * step. Output is only available at time_end. 
   *
   * @param grid a reference to a Grid object
   * @param time_step current time step
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  void calculate_result(const Grid &grid, 
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

#endif // POINT_DFT_RESULT_H
