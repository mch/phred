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

#ifndef UPML_BC_H
#define UPML_BC_H

#include <iostream>
#include <fstream>
using namespace std;

#include "BoundaryCondition.hh"
#include "UPmlCommon.hh"

/**
 * Uniaxial PML implementation, per Gedney. 
 *
 * Implmentation based on chapter 5 of "Advances in Engineering
 * Electromagnetics," Allen Taflove ed., Artech House, 1998
 */ 
class UPml : public BoundaryCond
{
private:
protected:
  UPmlCommon *common_;

  /**
   * Stored flux values from the last time step, used to compute a
   * second-order accurate frequency dependence.  
   */ 
  field_t *dx_, *dy_, *dz_;
  field_t *bx_, *by_, *bz_;

  /**
   * Additional auxiliary variables used by other dispersions
   */ 
  field_t *aux1_x_, *aux1_y_, *aux1_z_;

  /**
   * Sigma max, if the user choses to set it. When set to zero (the
   * default), sigma max = sigma opt = (m+1)/(150pi sqrt(eps_r) dx)
   */
  mat_coef_t sigma_max_;

  /**
   * Polynomial order for calculating the conductivity profile of the
   * UPML. Defaults to 4. 
   */ 
  unsigned int poly_order_;

  /**
   * The reletive permittivity to be used for the purposes of
   * calculating the maximum conductivity in the PML. If there are big
   * permittivity continuities withing the PML, the best performace
   * may be to average them. Defaults to 1.0. 
   */ 
  mat_coef_t eps_opt_;

  /**
   * The ratio between sigma_opt, calculated from the polynomial order
   * and material properties, to sigma_max, the maximum condictivity
   * of the UPML. Defaults to 1.0. 
   */ 
  mat_coef_t sigma_ratio_;

  // Update equations
  void update_ex(Grid &grid);
  void update_ey(Grid &grid);
  void update_ez(Grid &grid);

  void update_hx(Grid &grid);
  void update_hy(Grid &grid);
  void update_hz(Grid &grid);

  void compute_regions(Face face, const Grid &grid);

public:
  UPml();
  ~UPml();

  /**
   * Applies a boundary condition to a face of the grid.
   *
   * @param face the face to apply the boundary condition to. 
   * @param grid the grid to apply the boundary condition to. 
   * @param the field components to affect. 
   */
  virtual void apply(Face face,
                     Grid &grid, 
                     FieldType type);

  /**
   * Boundary condition type. Subclasses may implement. Only really
   * important for PMLs.
   */
  virtual BoundaryCondition get_type() const
  {
    return UPML;
  }

  /**
   * Init the boundary condition. 
   */ 
  virtual void init(const Grid &grid, Face face);

  /**
   * Cleanup the boundary condition. 
   */ 
  virtual void deinit(const Grid &grid, Face face);

  /**
   * Transfer 
   * 
   * @param sd the subdomain boundary condition that has to exchange
   * the data. 
   * @param bcface the face the PML is on
   * @param sdface the face the subdmoain is on
   */ 
  virtual void add_sd_bcs(SubdomainBc *sd, Face bcface, Face sdface);

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

  inline mat_coef_t get_sigma_ratio()
  { return sigma_ratio_; } 

  inline void set_sigma_ratio(mat_coef_t sr)
  { 
    if (sr > 0.0)
      sigma_ratio_ = sr; 
  }
};

#endif // UPML_BC_H
