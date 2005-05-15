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


// ADE Drude was the first one implemented, but the interior grid uses
// a drude forumulation that was developed using the Z transform. Use
// of the ADE in the UPML and the Z transform in the interior seems to
// cause instability (or there is an error in one of the
// constants). The Z transform formulation should also be more
// accurate.
#undef ADE_DRUDE
// #define ADE_DRUDE 1


/**
 * A class that contains data that is needed for all Uniaxial PML's on
 * a grid. This is a singleton; the assumption is that there is only
 * ever one grid being used in a program at a time.
 *
 * This class is cabable of matching interiors with inhomogenous
 * mixtures of dielectrics, lossy dielectrics, Drude (non magnetic
 * plasma) material. 
 *
 * Lorentz and Debye materials coming soon. 
 */ 
class UPmlCommon
{
private:
  UPmlCommon(const Grid &grid);
protected:
  const Grid &grid_;

  bool inited_; /**< Set to true once init_coeffs has been
                   called. Prevents multiple init_coeffs calls. */ 

  unsigned int poly_order_; /**< Order of the polynomial used to shape
                               the conductivity */

  mat_coef_t sigma_max_; /**< Sigma max; use 0 to calculate */

  /**
   * The reletive permittivity to be used for the purposes of
   * calculating the maximum conductivity in the PML. If there are big
   * permittivity continuities withing the PML, the best performace
   * may be to average them.
   */ 
  mat_coef_t eps_opt_;

  /**
   * The ratio between sigma_opt, calculated from the polynomial order
   * and material properties, to sigma_max, the maximum condictivity
   * of the UPML. Defaults to 1.0. 
   */ 
  mat_coef_t sigma_ratio_;

  /**
   * This parameter can be set to a value > 1 to help terminate
   * evenescant waves. It is graded into the PML like sigma max. 
   */ 
  mat_coef_t k_max_;

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

  // Material type, i.e. the dispersion to invoke
  MaterialType *mtype_;

  // Constants for the lossy dispersion, one for each material. 
  float *lossyA_, *lossyB_;

  // Constants for the Drude (unmagnitized plasma) dispersion, one for
  // each material, so they can be indexed using mid. 
#ifdef ADE_DRUDE
  float *drudeC1_, *drudeC2_, *drudeC3_, *drudeC4_, *drudeC5_;
#else
  float *vcdt_, *omegasq_;
#endif

  /**
   * Debye material constants, one for each material
   */ 
  field_t *debyeA_, *debyeB_, *debyeC_;

  void init_constants();

  void free_sigmas();
  void init_sigmas();

  mat_coef_t calc_sigma_max(delta_t delta);

public:
  ~UPmlCommon();

  static UPmlCommon *get_upml_common(Grid &grid);

  void init_coeffs();

  /**
   * Called by UPml::deinit. Free memory and reset the object.
   */ 
  void deinit();

  inline float Ax(loop_idx_t i) const
  {
    assert(i < x_size_);
    return Ax_[i];
  }

  inline float Ay(loop_idx_t i) const
  {
    assert(i < y_size_);
    return Ay_[i];
  }

  inline float Az(loop_idx_t i) const
  {
    assert(i < z_size_);
    return Az_[i];
  }

  inline float Bx(loop_idx_t i) const
  {
    assert(i < x_size_);
    return Bx_[i];
  }

  inline float By(loop_idx_t i) const
  {
    assert(i < y_size_);
    return By_[i];
  }

  inline float Bz(loop_idx_t i) const
  {
    assert(i < z_size_);
    return Bz_[i];
  }

  inline float Cx(loop_idx_t i) const
  {
    assert(i < x_size_);
    return Cx_[i];
  }

  inline float Cy(loop_idx_t i) const
  {
    assert(i < y_size_);
    return Cy_[i];
  }

  inline float Cz(loop_idx_t i) const
  {
    assert(i < z_size_);
    return Cz_[i];
  }

  inline float Dx(loop_idx_t i) const
  {
    assert(i < x_size_);
    return Dx_[i];
  }

  inline float Dy(loop_idx_t i) const
  {
    assert(i < y_size_);
    return Dy_[i];
  }

  inline float Dz(loop_idx_t i) const
  {
    assert(i < z_size_);
    return Dz_[i];
  }

  inline float er(mat_idx_t mid) const
  {
    return er_[mid];
  }

  inline float ur(mat_idx_t mid) const
  {
    return ur_[mid];
  }

  inline float lossy_A(mat_idx_t mid) const
  { return lossyA_[mid]; }

  inline float lossy_B(mat_idx_t mid) const
  { return lossyB_[mid]; }

  inline MaterialType mtype(mat_idx_t mid) const
  { return mtype_[mid]; }

  // ADE Drude model constants
#ifdef ADE_DRUDE
  inline float drude_c1(mat_idx_t mid) const
  { return drudeC1_[mid]; }

  inline float drude_c2(mat_idx_t mid) const
  { return drudeC2_[mid]; }

  inline float drude_c3(mat_idx_t mid) const
  { return drudeC3_[mid]; }

  inline float drude_c4(mat_idx_t mid) const
  { return drudeC4_[mid]; }

  inline float drude_c5(mat_idx_t mid) const
  { return drudeC5_[mid]; }

#else
  // Z transform Drude model constants
  inline float get_vcdt(mat_idx_t mid) const
  { return vcdt_[mid]; }

  inline float get_omegasq(mat_idx_t mid) const
  { return omegasq_[mid]; }
#endif

  // Debye constants
  inline field_t debyeA(mat_idx_t mid) const
  { return debyeA_[mid]; }

  inline field_t debyeB(mat_idx_t mid) const
  { return debyeB_[mid]; }

  inline field_t debyeC(mat_idx_t mid) const
  { return debyeC_[mid]; }

  // Parameters set by the UPml object set by the user
  inline unsigned int get_poly_order() const
  { return poly_order_; }

  inline void set_poly_order(unsigned int p) 
  { poly_order_ = p; }

  inline mat_coef_t get_sigma_max() const
  { return sigma_max_; }

  inline void set_sigma_max(mat_coef_t sm) 
  { sigma_max_ = sm; }

  inline mat_coef_t get_eps_opt() const
  { return eps_opt_; }

  inline void set_eps_opt(mat_coef_t eps_opt) 
  { eps_opt_ = eps_opt; }

  inline mat_coef_t set_sigma_ratio()
  { return sigma_ratio_; } 

  inline void set_sigma_ratio(mat_coef_t sr)
  { 
    if (sr > 0.0)
      sigma_ratio_ = sr; 
  }

  inline void set_k_max(mat_coef_t km)
  { k_max_ = km; }

  inline mat_coef_t get_k_max()
  { return k_max_; }

};

#endif // UPML_COMMON_H
