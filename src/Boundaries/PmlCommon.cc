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

#include "PmlCommon.hh"
#include "Grid.hh"
#include "Exceptions.hh"

#include <math.h>

PmlCommon::PmlCommon(const Grid &grid)
  : grid_(grid), ratio_x_(0), ratio_star_x_(0), 
    ratio_y_(0), ratio_star_y_(0), 
    ratio_z_(0), ratio_star_z_(0), 
    e_x_coef1_(0), e_x_coef2_(0), 
    e_y_coef1_(0), e_y_coef2_(0), 
    e_z_coef1_(0), e_z_coef2_(0), 
    h_x_coef1_(0), h_x_coef2_(0), 
    h_y_coef1_(0), h_y_coef2_(0), 
    h_z_coef1_(0), h_z_coef2_(0)
{}

PmlCommon::~PmlCommon()
{
  free_coeffs();
}

PmlCommon *PmlCommon::get_pml_common(Grid &grid)
{
  if (grid.get_define_mode())
    return 0;

  void *tmp = grid.get_auxdata(PML_COMMON);

  if (tmp)
    return (PmlCommon *)tmp;

  PmlCommon *pml_common = new PmlCommon(grid);
  if (!pml_common)
    return 0; // stack unwinding may be too costly here

  pml_common->init_coeffs(grid);

  return pml_common;
}

void PmlCommon::alloc_coeffs(const Grid &grid)
{
  free_coeffs();

  ratio_x_ = new float[grid.get_ldx_sd()];
  ratio_star_x_ = new float[grid.get_ldx_sd()];
  
  ratio_y_ = new float[grid.get_ldy_sd()];
  ratio_star_y_ = new float[grid.get_ldy_sd()];
  
  ratio_z_ = new float[grid.get_ldz_sd()];
  ratio_star_z_ = new float[grid.get_ldz_sd()];
  
  e_x_coef1_ = new float[grid.get_ldx_sd()];
  e_x_coef2_ = new float[grid.get_ldx_sd()];
  
  e_y_coef1_ = new float[grid.get_ldy_sd()];
  e_y_coef2_ = new float[grid.get_ldy_sd()];
  
  e_z_coef1_ = new float[grid.get_ldz_sd()];
  e_z_coef2_ = new float[grid.get_ldz_sd()];
  
  h_x_coef1_ = new float[grid.get_ldx_sd()];
  h_x_coef2_ = new float[grid.get_ldx_sd()];
  
  h_y_coef1_ = new float[grid.get_ldy_sd()];
  h_y_coef2_ = new float[grid.get_ldy_sd()];
  
  h_z_coef1_ = new float[grid.get_ldz_sd()];
  h_z_coef2_ = new float[grid.get_ldz_sd()];  

  if (!ratio_x_ || !ratio_star_x_ || !ratio_y_ || !ratio_star_y_ || 
      !ratio_z_ || !ratio_star_z_ || !e_x_coef1_ || !e_x_coef2_ ||
      !e_y_coef1_ || !e_y_coef2_ || !e_z_coef1_ || !e_z_coef2_ ||
      !h_x_coef1_ || !h_x_coef2_ || !h_y_coef1_ || !h_y_coef2_ || 
      !h_z_coef1_ || !h_z_coef2_)
  {
    free_coeffs();
    throw MemoryException();
  }

  memset(ratio_x_, 0, sizeof(float) * grid.get_ldx_sd());
  memset(ratio_star_x_, 0, sizeof(float) * grid.get_ldx_sd());

  memset(ratio_y_, 0, sizeof(float) * grid.get_ldy_sd());
  memset(ratio_star_y_, 0, sizeof(float) * grid.get_ldy_sd());

  memset(ratio_z_, 0, sizeof(float) * grid.get_ldz_sd());
  memset(ratio_star_z_, 0, sizeof(float) * grid.get_ldz_sd());

  memset(e_x_coef1_, 0, sizeof(float) * grid.get_ldx_sd());
  memset(e_x_coef2_, 0, sizeof(float) * grid.get_ldx_sd());

  memset(e_y_coef1_, 0, sizeof(float) * grid.get_ldy_sd());
  memset(e_y_coef2_, 0, sizeof(float) * grid.get_ldy_sd());

  memset(e_z_coef1_, 0, sizeof(float) * grid.get_ldz_sd());
  memset(e_z_coef2_, 0, sizeof(float) * grid.get_ldz_sd());

  memset(h_x_coef1_, 0, sizeof(float) * grid.get_ldx_sd());
  memset(h_x_coef2_, 0, sizeof(float) * grid.get_ldx_sd());

  memset(h_y_coef1_, 0, sizeof(float) * grid.get_ldy_sd());
  memset(h_y_coef2_, 0, sizeof(float) * grid.get_ldy_sd());

  memset(h_z_coef1_, 0, sizeof(float) * grid.get_ldz_sd());
  memset(h_z_coef2_, 0, sizeof(float) * grid.get_ldz_sd());

#ifndef NDEBUG
  dimx_ = grid.get_ldx_sd();
  dimy_ = grid.get_ldy_sd();
  dimz_ = grid.get_ldz_sd();
  cerr << "Set dimx_ to " << dimx_ << ", dimy_ to " << dimy_
       << ", and dimz_ to " << dimz_ << endl;
#endif
}

void PmlCommon::free_coeffs()
{
  if (ratio_x_)
    delete[] ratio_x_;

  if (ratio_star_x_)
    delete[] ratio_star_x_;
  
  if (ratio_y_)
    delete[] ratio_y_;

  if (ratio_star_y_)
    delete[] ratio_star_y_;
  
  if (ratio_z_)
    delete[] ratio_z_;
  if (ratio_star_z_)
    delete[] ratio_star_z_;
  
  if (e_x_coef1_)
    delete[] e_x_coef1_;
  if (e_x_coef2_)
    delete[] e_x_coef2_;
  
  if (e_y_coef1_)
    delete[] e_y_coef1_;
  if (e_y_coef2_)
    delete[] e_y_coef2_;
  
  if (e_z_coef1_)
    delete[] e_z_coef1_;
  if (e_z_coef2_)
    delete[] e_z_coef2_;
  
  if (h_x_coef1_)
    delete[] h_x_coef1_;
  if (h_x_coef2_)
    delete[] h_x_coef2_;
  
  if (h_y_coef1_)
    delete[] h_y_coef1_;
  if (h_y_coef2_)
    delete[] h_y_coef2_;
  
  if (h_z_coef1_)
    delete[] h_z_coef1_;
  if (h_z_coef2_)
    delete[] h_z_coef2_;

  ratio_x_ = ratio_star_x_ = ratio_y_ = ratio_star_y_ = 0;
  ratio_z_ = ratio_star_z_ = 0;
  e_x_coef1_ = e_x_coef2_ = e_y_coef1_ = e_y_coef2_ = 0;
  e_z_coef1_ = e_z_coef2_ = h_x_coef1_ = h_x_coef2_ = 0;
  h_y_coef1_ = h_y_coef2_ = h_z_coef1_ = h_z_coef2_ = 0;
}


void PmlCommon::init_coeffs(Grid &grid) 
{
  alloc_coeffs(grid);

  GridInfo &gi = grid.get_grid_info();

  for (int i = 0; i < 6; i++)
  {
    if (gi.get_bc_type(static_cast<Face>(i)) == PML) {
      BoundaryCond &bc = gi.get_boundary(static_cast<Face>(i));
      Pml *p = dynamic_cast<Pml *>(&bc);

      if (p) 
      {
        init_ratios(static_cast<Face>(i), grid, p);
      } // if (p)
    }
  } // for


//   cout << endl << "item\tratio_x\t\tratio_star_x " << endl;
//   for (int j = 0; j < grid.get_ldx(); j++)
//     cout << j << "\t" << ratio_x_[j] << "\t\t" << ratio_star_x_[j] << endl;

//   cout << endl << "item\tratio_y\t\tratio_star_y " << endl;
//   for (int j = 0; j < grid.get_ldy(); j++)
//     cout << j << "\t" << ratio_y_[j] << "\t\t" << ratio_star_y_[j] << endl;

//   cout << endl << "item\tratio_z\t\tratio_star_z " << endl;
//   for (int j = 0; j < grid.get_ldz(); j++)
//     cout << j << "\t" << ratio_z_[j] << "\t\t" << ratio_star_z_[j] << endl;
  

  // And finally, init coefficients
  delta_t dt = grid.get_deltat();

  for(unsigned int i = 0; i < grid.get_ldx_sd(); i++)
  {
    e_x_coef1_[i] = (ratio_x_[i] == 0.0 ? 1.0
                     : exp(-ratio_x_[i] * dt));

    e_x_coef2_[i] = (ratio_x_[i] == 0.0 ? 1.0
                     : (1.0 - exp(-ratio_x_[i] * dt))
                     / (ratio_x_[i] * dt));
                        
    h_x_coef1_[i] = (ratio_star_x_[i] == 0.0 ? 1.0 
                     : exp(-ratio_star_x_[i] * dt));

    h_x_coef2_[i] = (ratio_star_x_[i] == 0.0 ? 1.0 
                     : (1.0 - exp(-ratio_star_x_[i] * dt))
                     / (ratio_star_x_[i] * dt));
  }
                

  for(unsigned int j = 0; j < grid.get_ldy_sd(); j++)
  {
    e_y_coef1_[j] = (ratio_y_[j] == 0.0 ? 1.0 
                     : exp(-ratio_y_[j] * dt));

    e_y_coef2_[j] = (ratio_y_[j] == 0.0 ? 1.0 
                     : (1.0 - exp(-ratio_y_[j] * dt))
                     / (ratio_y_[j] * dt));
                        
    h_y_coef1_[j] = (ratio_star_y_[j] == 0.0 ? 1.0 
                     : exp(-ratio_star_y_[j] * dt));

    h_y_coef2_[j] = (ratio_star_y_[j] == 0.0 ? 1.0 
                     : (1.0 - exp(-ratio_star_y_[j] * dt))
                     / (ratio_star_y_[j] * dt));
  }

  for(unsigned int k = 0; k < grid.get_ldz_sd(); k++)
  {
    e_z_coef1_[k] = (ratio_z_[k] == 0.0 ? 1.0
                     : exp(-ratio_z_[k] * dt));

    e_z_coef2_[k] = (ratio_z_[k] == 0.0 ? 1. 
                     : (1.0 - exp(-ratio_z_[k] * dt))
                     / (ratio_z_[k] * dt));
                        
    h_z_coef1_[k] = (ratio_star_z_[k] == 0.0 ? 1.0 
                     : exp(-ratio_star_z_[k] * dt));

    h_z_coef2_[k] = (ratio_star_z_[k] == 0.0 ? 1.0 
                     : (1.0 - exp(-ratio_star_z_[k] * dt))
                     / (ratio_star_z_[k] * dt));
  }
                  
//   cout << endl << "item\te_x_coef1\t\te_x_coef2\t\th_x_coef1\t\th_x_coef2 " << endl;
//   for (int j = 0; j < grid.get_ldx_sd(); j++)
//     cout << j << "\t" << e_x_coef1_[j] << "\t\t" << e_x_coef2_[j]
//          << "\t\t" << h_x_coef1_[j] << "\t\t" << h_x_coef2_[j] << endl;

//   cout << endl << "item\te_y_coef1\t\te_y_coef2\t\th_y_coef1\t\th_y_coef2 " << endl;
//   for (int j = 0; j < grid.get_ldy_sd(); j++)
//     cout << j << "\t" << e_y_coef1_[j] << "\t\t" << e_y_coef2_[j]
//          << "\t\t" << h_y_coef1_[j] << "\t\t" << h_y_coef2_[j] << endl;

//   cout << endl << "item\te_z_coef1\t\te_z_coef2\t\th_z_coef1\t\th_z_coef2 " << endl;
//   for (int j = 0; j < grid.get_ldz_sd(); j++)
//     cout << j << "\t" << e_z_coef1_[j] << "\t\t" << e_z_coef2_[j]
//          << "\t\t" << h_z_coef1_[j] << "\t\t" << h_z_coef2_[j] << endl;
  
}

void PmlCommon::init_ratios(Face face, const Grid &grid, Pml *p)
{
  delta_t d_space;
  int step; 
  unsigned int idx1, idx2;
  float *arr1, *arr2;

  switch (face)
  {
  case BACK:
    d_space = grid.get_deltax();
    step = 1;
    idx1 = 0;
    idx2 = 0;
    arr1 = ratio_x_;
    arr2 = ratio_star_x_;
    break;

  case FRONT:
    d_space = grid.get_deltax();
    step = -1; 
    idx1 = grid.get_ldx_sd() - 1;
    idx2 = grid.get_ldx_sd() - 2;
    arr1 = ratio_x_;
    arr2 = ratio_star_x_;
    break;

  case LEFT:
    d_space = grid.get_deltay();
    step = 1;
    idx1 = 0;
    idx2 = 0;
    arr1 = ratio_y_;
    arr2 = ratio_star_y_;
    break;

  case RIGHT:
    d_space = grid.get_deltay();
    step = -1; 
    idx1 = grid.get_ldy_sd() - 1;
    idx2 = grid.get_ldy_sd() - 2;
    arr1 = ratio_y_;
    arr2 = ratio_star_y_;
    break;

  case BOTTOM:
    d_space = grid.get_deltaz();
    step = 1;
    idx1 = 0;
    idx2 = 0;
    arr1 = ratio_z_;
    arr2 = ratio_star_z_;
    break;

  case TOP:
    d_space = grid.get_deltaz();
    step = -1; 
    idx1 = grid.get_ldz_sd() - 1;
    idx2 = grid.get_ldz_sd() - 2;
    arr1 = ratio_z_;
    arr2 = ratio_star_z_;
    break;
  }

  for (unsigned int i = 0; i <= p->get_thickness(); i++)
  {
    arr1[idx1] = 1.0 / d_space 
      * ( p->sigma_over_eps_int((i+0.5)*d_space) 
          - p->sigma_over_eps_int((i-0.5)*d_space));
    
    arr2[idx2] = 1.0 / d_space 
      * ( p->sigma_over_eps_int((i+1.0)*d_space) 
          - p->sigma_over_eps_int((i)*d_space));
    
    idx1 += step;
    idx2 += step;
  }

}

