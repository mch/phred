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

#ifndef PML_COMMON_H
#define PML_COMMON_H

#include "Types.hh"
#include "config.h"

#include <assert.h>

class Pml;
class Grid;

/**
 * Data common to a set of PML's. This holds coefficients and stuff
 * that are used by all PML's. This is intended to be a member of a
 * Grid class, which will call the alloc_coeffs and free_coeffs
 * methods. Clients other the a Grid or a Pml won't be able to do much
 * with this...
 */
class PmlCommon
{
  //friend class Pml;
  //friend class Grid;

private:
  const Grid &grid_;

  float *ratio_x_;
  float *ratio_star_x_;
  
  float *ratio_y_;
  float *ratio_star_y_;
  
  float *ratio_z_;
  float *ratio_star_z_;
  
  float *e_x_coef1_;
  float *e_x_coef2_;
  
  float *e_y_coef1_;
  float *e_y_coef2_;
  
  float *e_z_coef1_;
  float *e_z_coef2_;
  
  float *h_x_coef1_;
  float *h_x_coef2_;
  
  float *h_y_coef1_;
  float *h_y_coef2_;
  
  float *h_z_coef1_;
  float *h_z_coef2_;

  // For assertions, so we know when we step out of bounds
  unsigned int dimx_;
  unsigned int dimy_;
  unsigned int dimz_;

  /**
   * Setup the coefficients 
   */
  void alloc_coeffs(const Grid &grid);

  /**
   * Free the coefficients.
   */
  void free_coeffs();

  /**
   * Helper function to initialize ratios
   */
  void init_ratios(Face face, const Grid &grid, Pml *p);

  /**
   * Constructor
   */
  PmlCommon(const Grid &grid);

  /**
   * Set the common PML parameters and calculate coeffs and
   * stuff. Called by Grid when leaving define mode, when all the
   * boundary conditions have been set. 
   *
   * @param grid The grid the PML is being applied to. 
   */
  void init_coeffs(Grid &grid);

public:
  /**
   * Destructornator!
   */
  ~PmlCommon();

  /**
   * Returns the only available instance of this class. If the grid is
   * in define mode, null will be returned.
   *
   * @param grid required to initialize the class if this is the first
   * time the instance has been requrested
   *
   * @return If the grid is in define mode, null will be returned,
   * otherwise, the pml common object will be returned. 
   */
  static PmlCommon *get_pml_common(Grid &grid);
  
  /**
   * Return a e_x_coef1 coefficient
   */
  inline float get_e_x_coef1(unsigned int i)
  {
    assert(i < dimx_);
    return e_x_coef1_[i];
  }

  /**
   * Return a e_y_coef1 coefficient
   */
  inline float get_e_y_coef1(unsigned int j)
  {
    assert(j < dimy_);
    return e_y_coef1_[j];
  }

  /**
   * Return a e_z_coef1 coefficient
   */
  inline float get_e_z_coef1(unsigned int k)
  {
    assert(k < dimz_);
    return e_z_coef1_[k];
  }

  /**
   * Return a e_x_coef2 coefficient
   */
  inline float get_e_x_coef2(unsigned int i)
  {
    assert(i < dimx_);
    return e_x_coef2_[i];
  }

  /**
   * Return a e_y_coef2 coefficient
   */
  inline float get_e_y_coef2(unsigned int j)
  {
    assert(j < dimy_);
    return e_y_coef2_[j];
  }

  /**
   * Return a e_z_coef2 coefficient
   */
  inline float get_e_z_coef2(unsigned int k)
  {
    assert(k < dimz_);
    return e_z_coef2_[k];
  }

  /**
   * Return a h_x_coef1 coefficient
   */
  inline float get_h_x_coef1(unsigned int i)
  {
    assert(i < dimx_);
    return h_x_coef1_[i];
  }

  /**
   * Return a h_y_coef1 coefficient
   */
  inline float get_h_y_coef1(unsigned int j)
  {
    assert(j < dimy_);
    return h_y_coef1_[j];
  }

  /**
   * Return a h_z_coef1 coefficient
   */
  inline float get_h_z_coef1(unsigned int k)
  {
    assert(k < dimz_);
    return h_z_coef1_[k];
  }

  /**
   * Return a h_x_coef2 coefficient
   */
  inline float get_h_x_coef2(unsigned int i)
  {
    assert(i < dimx_);
    return h_x_coef2_[i];
  }

  /**
   * Return a h_y_coef2 coefficient
   */
  inline float get_h_y_coef2(unsigned int j)
  {
    assert(j < dimy_);
    return h_y_coef2_[j];
  }

  /**
   * Return a h_z_coef2 coefficient
   */
  inline float get_h_z_coef2(unsigned int k)
  {
    assert(k < dimz_);
    return h_z_coef2_[k];
  }

};

#endif
