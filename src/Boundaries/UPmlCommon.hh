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

#ifndef UPML_COMMON_H
#define UPML_COMMON_H

#include "Grid.hh"

/**
 * A class that contains data that is needed for all Uniaxial PML's on
 * a grid. This is a singleton; the assumption is that there is only
 * ever one grid being used in a program at a time.
 */ 
class UPmlCommon
{
private:
  UPmlCommon(const Grid &grid);
protected:
  const Grid &grid_;

  unsigned int poly_order_; /**< Order of the polynomial used to shape
                               the conductivity */

  unsigned int thicknesses_[6]; /**< UPml thicknesses */

  /**
   * Conductivities. 
   */
  float *sigma_x_;
  float *sigma_y_;
  float *sigma_z_;

  /**
   * 4d coefficient matricies. Only access these using the coeff_id
   * function. 
   */ 
  float *c1_;
  float *c2_;
  float *c3_;
  float *c4_;
  float *c5_;
  float *c6_;

  unsigned int num_materials_;
  unsigned int x_size_; /**< Length in X of coefficient matricies */
  unsigned int y_size_;
  unsigned int z_size_;

  // The above sizes are the sum of the thicknesses of the pml at
  // opposite faces plus one for the regions where there is no
  // overlap. 

  /**
   * Returns a sigma value along the x axis.
   *
   * @param mid material id
   * @param x x coordinate
   */
  inline float sigma_x(unsigned int mid, unsigned int x)
  {
    unsigned int xcoord = 0; 
    float ret = 0.0;

    if (x < thicknesses_[BACK])
      xcoord = x;
    if (x > (grid_.get_ldx() - thicknesses_[FRONT]))
      xcoord = thicknesses_[BACK]
        + (x - (grid_.get_ldx() - thicknesses_[FRONT]));

    if (mid > 0)
      ret = sigma_x_[xcoord + ((mid - 1) * num_materials_)];

    return ret;
  }

  /**
   * Returns a sigma value along the y axis.
   *
   * @param mid material id
   * @param y y coordinate
   */
  inline float sigma_y(unsigned int mid, unsigned int y)
  {
    unsigned int ycoord = 0; 
    float ret = 0.0;

    if (y < thicknesses_[LEFT])
      ycoord = y;
    if (y > (grid_.get_ldy() - thicknesses_[RIGHT]))
      ycoord = thicknesses_[LEFT]
        + (y - (grid_.get_ldy() - thicknesses_[RIGHT]));

    if (mid > 0)
      ret = sigma_y_[ycoord + ((mid - 1) * num_materials_)];

    return ret;
  }

  /**
   * Returns a sigma value along the z axis.
   *
   * @param mid material id
   * @param z z coordinate
   */
  inline float sigma_z(unsigned int mid, unsigned int z)
  {
    unsigned int zcoord = 0; 
    float ret = 0.0;

    if (z < thicknesses_[BOTTOM])
      zcoord = z;
    if (z > (grid_.get_ldz() - thicknesses_[TOP]))
      zcoord = thicknesses_[BOTTOM]
        + (z - (grid_.get_ldz() - thicknesses_[TOP]));

    if (mid > 0)
      ret = sigma_z_[zcoord + ((mid - 1) * num_materials_)];

    return ret;
  }

  /**
   * This function computes the index required to retrieve the
   * correct coefficient. 
   */ 
  inline unsigned int coeff_id(unsigned int material_id, 
                               unsigned int x, unsigned int y, 
                               unsigned int z)
  {
    return 0; 
  }

  void init_coeffs(Grid &grid);

  void free_sigmas();
  void init_sigmas();

  mat_coef_t calc_sigma_max(mat_prop_t eps, delta_t delta);

public:
  ~UPmlCommon();

  static UPmlCommon *get_upml_common(Grid &grid);
};

#endif // UPML_COMMON_H
