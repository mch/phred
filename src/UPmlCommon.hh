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

/**
 * A class that contains data that is needed for all Uniaxial PML's on
 * a grid. This is a singleton; the assumption is that there is only
 * ever one grid being used in a program at a time.
 */ 
class UPmlCommon
{
private:
  UPmlCommon();
protected:

  /**
   * Conductivities. 
   */
  float **sigma_x_;
  float **sigma_y_;
  float **sigma_z_;

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
   * This function computes the index into the sigma arrays
   */
  inline unsigned int sigma_id();

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
  void init_sigmas(Face face, const Grid &grid, UPml *p);

public:
  ~UPmlCommon();

  static UPmlCommon *get_upml_common(Grid &grid);
};

#endif // UPML_COMMON_H
