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

#ifndef UPML_COMMON_H
#define UPML_COMMON_H

#include "../Grid.hh"

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

  // Quick and dirty, these will be equal to the local grid
  // size. Since a lot of the interior space will be equal to zero,
  // this should be changed so that it's 2*PML thickness.
  unsigned int x_size_; /**< Length in X of coefficient matricies */
  unsigned int y_size_; /**< Length in Y of coefficient matricies */
  unsigned int z_size_; /**< Length in Z of coefficient matricies */

  /**
   * Conductivity along the X axis inside the PML. 
   */
  float *sigma_x_;

  /**
   * Conductivity along the Y axis inside the PML. 
   */
  float *sigma_y_;

  /**
   * Conductivity along the Z axis inside the PML. 
   */
  float *sigma_z_;

  // Evenescant wave attenuation term
  float *kx_, *ky_, *kz_;

  // Ax(x) = ( 2*eps_0*Kx(x) - dt*sigma_x(x) ) 
  //         / ( 2*eps_0*Kx(x) + dt*sigma_x(x) )
  // etc.
  float *Ax_, *Ay_, *Az_;

  // Bx(x) = 2*eps_0*dt
  //         / ( 2*eps_0*Kx(x) + dt*sigma_x(x) )
  //         * (dx / (dy*dz)) ???????? CHECK
  // etc.
  float *Bx_, *By_, *Bz_;

  // Cx(x) = ( 2*eps_0*Kx(x) + dt*sigma_x(x) ) 
  //         / 2*eps_0*dt
  float *Cx_, *Cy_, *Cz_;

  // Dx(x) = ( 2*eps_0*Kx(x) - dt*sigma_x(x) )
  //         / 2*eps_0*dt
  float *Dx_, *Dy_, *Dz_;

  // Inverse permittivity, one for each material
  float *er_;

  // Inverse permeability, one for each material
  float *ur_;

//   /**
//    * Returns a sigma value along the x axis.
//    *
//    * @param mid material id
//    * @param x x coordinate
//    */
//   inline float sigma_x(unsigned int mid, unsigned int x)
//   {
//     unsigned int xcoord = 0; 
//     float ret = 0.0;

//     if (x < thicknesses_[BACK])
//       xcoord = x;
//     if (x > (grid_.get_ldx() - thicknesses_[FRONT]))
//       xcoord = thicknesses_[BACK]
//         + (x - (grid_.get_ldx() - thicknesses_[FRONT]));

//     if (mid > 0)
//       ret = sigma_x_[xcoord + ((mid - 1) * num_materials_)];

//     return ret;
//   }

//   /**
//    * Returns a sigma value along the y axis.
//    *
//    * @param mid material id
//    * @param y y coordinate
//    */
//   inline float sigma_y(unsigned int mid, unsigned int y)
//   {
//     unsigned int ycoord = 0; 
//     float ret = 0.0;

//     if (y < thicknesses_[LEFT])
//       ycoord = y;
//     if (y > (grid_.get_ldy() - thicknesses_[RIGHT]))
//       ycoord = thicknesses_[LEFT]
//         + (y - (grid_.get_ldy() - thicknesses_[RIGHT]));

//     if (mid > 0)
//       ret = sigma_y_[ycoord + ((mid - 1) * num_materials_)];

//     return ret;
//   }

//   /**
//    * Returns a sigma value along the z axis.
//    *
//    * @param mid material id
//    * @param z z coordinate
//    */
//   inline float sigma_z(unsigned int mid, unsigned int z)
//   {
//     unsigned int zcoord = 0; 
//     float ret = 0.0;

//     if (z < thicknesses_[BOTTOM])
//       zcoord = z;
//     if (z > (grid_.get_ldz() - thicknesses_[TOP]))
//       zcoord = thicknesses_[BOTTOM]
//         + (z - (grid_.get_ldz() - thicknesses_[TOP]));

//     if (mid > 0)
//       ret = sigma_z_[zcoord + ((mid - 1) * num_materials_)];

//     return ret;
//   }

  void init_coeffs();
  void init_constants();

  void free_sigmas();
  void init_sigmas();

  mat_coef_t calc_sigma_max(mat_prop_t eps, delta_t delta);

public:
  ~UPmlCommon();

  static UPmlCommon *get_upml_common(Grid &grid);

  inline const float Ax(loop_idx_t i)
  {
    return Ax_[i];
  }

  inline const float Ay(loop_idx_t i)
  {
    return Ay_[i];
  }

  inline const float Az(loop_idx_t i)
  {
    return Az_[i];
  }

  inline const float Bx(loop_idx_t i)
  {
    return Bx_[i];
  }

  inline const float By(loop_idx_t i)
  {
    return By_[i];
  }

  inline const float Bz(loop_idx_t i)
  {
    return Bz_[i];
  }

  inline const float Cx(loop_idx_t i)
  {
    return Cx_[i];
  }

  inline const float Cy(loop_idx_t i)
  {
    return Cy_[i];
  }

  inline const float Cz(loop_idx_t i)
  {
    return Cz_[i];
  }

  inline const float Dx(loop_idx_t i)
  {
    return Dx_[i];
  }

  inline const float Dy(loop_idx_t i)
  {
    return Dy_[i];
  }

  inline const float Dz(loop_idx_t i)
  {
    return Dz_[i];
  }

  inline const float er(mat_idx_t mid)
  {
    return er_[mid];
  }

  inline const float ur(mat_idx_t mid)
  {
    return ur_[mid];
  }
};

#endif // UPML_COMMON_H
