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

#include "UPml.hh"
#include "Exceptions.hh"
#include "Grid.hh"
#include "Constants.hh"

#ifdef USE_OPENMP
#include <omp.h>
#endif

UPml::UPml()
  : d_(0), h_(0), sigmas_(0), poly_order_(4)
{

}

UPml::~UPml()
{

}

void UPml::init(const Grid &grid, Face face)
{
  compute_regions(face, grid);
  out_face_ = face;

  // Allocate memory
  alloc_coefs(grid.get_material_lib().num_materials());

  unsigned int sz = (grid_r_.xmax - grid_r_.xmin) 
    * (grid_r_.ymax - grid_r_.ymin) * (grid_r_.zmax - grid_r_.zmin);
  
  d_ = new field_t[sz];
  h_ = new field_t[sz];

  if (!d_ || !h_)
  {
    deinit(grid, face);
    throw MemoryException();
  }

  memset(d_, 0, sizeof(field_t) * sz);
  memset(h_, 0, sizeof(field_t) * sz);

  // Compute the material constants...
  field_t delta = 0;
  switch (face)
  {
  case FRONT:
  case BACK:
    delta = grid.get_deltax();
    break;
  case TOP:
  case BOTTOM:
    delta = grid.get_deltaz();
    break;
  case LEFT:
  case RIGHT:
    delta = grid.get_deltay();
    break;
  }

  mat_coef_t sigma_max = (poly_order_ + 1) 
    / (150 * PI * delta * sqrt(EPS_0));

  sigmas_ = new mat_coef_t[thickness_];

  for (int idx = 1; idx =< thickness_; idx++)
    sigmas_[idx] = sigma_max * pow(static_cast<float>(idx), 
                                   static_cast<float>(poly_order_)) 
      / pow(static_cast<float>(thickness_),
            static_cast<float>(poly_order_));

  // Build the shadow library. If the material's conductivity is
  // higher than the PML conductivity, leave it alone.
  MaterialLib *mlib = &(grid.get_material_lib());
  
  for (int idx = 0; idx < thicnkess_; idx++)
  {
    vector<Material>::iterator iter = mlib->get_material_iter_begin();
    vector<Material>::iterator iter_e = mlib->get_material_iter_end();
    int index = 0;

    Ca_[idx][index] = 1;
    Cbx_[idx][index] = Cby_[idx][index] = Cbz_[idx][index] = 0;

    Da_[index] = 1;
    Dbx_[idx][index] = Dby_[idx][index] = Dbz_[idx][index] = 0;
  
    ++index;
    while (iter != iter_e) 
    {
      mat_prop_t eps = (*iter).get_epsilon() * EPS_0;
      mat_prop_t sig = (*iter).get_sigma();
      mat_prop_t mu = (*iter).get_mu() * MU_0;
      mat_prop_t sigs = (*iter).get_sigma_star();

      if (sigmas_[idx] > sig)
        sig = sigmas_[idx];

      if (sigma_stars_[idx] > sigs)
        sigs = sigma_stars_[idx];

      Ca_[idx][index] = (1 - (sig * grid.get_deltat() * 0.5)/eps) / 
        (1 + (sig * grid.get_deltat() * 0.5)/eps);

      Da_[idx][index] = (1 - (sigs * grid.get_deltat() * 0.5)/mu) / 
        (1 + (sigs * grid.get_deltat() * 0.5)/mu);

    
      Cbx_[idx][index] = (grid.get_deltat() / (eps * grid.get_deltax())) / 
        (1 + (sig * grid.get_deltat() * 0.5)/eps);

      Dbx_[idx][index] = (grid.get_deltat() / (mu * grid.get_deltax())) / 
        (1 + (sigs * grid.get_deltat() * 0.5)/mu);

      Cby_[idx][index] = (grid.get_deltat() / (eps * grid.get_deltay())) / 
        (1 + (sig * grid.get_deltat() * 0.5)/eps);

      Dby_[idx][index] = (grid.get_deltat() / (mu * grid.get_deltay())) / 
        (1 + (sigs * grid.get_deltat() * 0.5)/mu);

      Cbz_[idx][index] = (grid.get_deltat() / (eps * grid.get_deltaz())) / 
        (1 + (sig * grid.get_deltat() * 0.5)/eps);

      Dbz_[idx][index] = (grid.get_deltat() / (mu * grid.get_deltaz())) / 
        (1 + (sigs * grid.get_deltat() * 0.5)/mu);

      C1_[idx][index] = (2 * eps * sigs - sig * grid.get_deltat())
        / (2 * eps * sigs + sig * grid.get_deltat());

      C2_[idx][index] = (2 * eps * grid.get_deltat()) 
        / (2 * eps * sigs + sig * grid.get_deltat());

      C3_[idx][index] = (2 * eps * 

      ++iter;
      ++index;      
    }
  }
}

void UPml::alloc_coefs(unsigned int num_materials)
{
  free_coefs();

  Ca_ = new mat_coef_t*[thickness_];
  Cbx_ = new mat_coef_t*[thickness_];
  Cby_ = new mat_coef_t*[thickness_];
  Cbz_ = new mat_coef_t*[thickness_];

  Da_ = new mat_coef_t*[thickness_];
  Dbx_ = new mat_coef_t*[thickness_];
  Dby_ = new mat_coef_t*[thickness_];
  Dbz_ = new mat_coef_t*[thickness_];
  
  C1_ = new mat_coef_t*[thicnkess_];
  C2_ = new mat_coef_t*[thicnkess_];
  C3_ = new mat_coef_t*[thicnkess_];
  C4_ = new mat_coef_t*[thicnkess_];
  C5_ = new mat_coef_t*[thicnkess_];
  C6_ = new mat_coef_t*[thicnkess_];

  if (!Ca_ || !Cbx_ || !Cbz_ || !Cby_ || !Da_ || !Dbx_ || !Dbz_ || !Dby_
      || !C1_ || !C2 || !C3 || !C4 || !C5 || !C6)
  {
    free_coefs();
    throw MemoryException();
  }

  for (int idx = 0; idx < thickness_; idx++)
  {
    Ca_[idx] = new mat_coef_t[num_materials];
    Cbx_[idx] = new mat_coef_t[num_materials];
    Cby_[idx] = new mat_coef_t[num_materials];
    Cbz_[idx] = new mat_coef_t[num_materials];

    Da_[idx] = new mat_coef_t[num_materials];
    Dbx_[idx] = new mat_coef_t[num_materials];
    Dby_[idx] = new mat_coef_t[num_materials];
    Dbz_[idx] = new mat_coef_t[num_materials];

    C1_[idx] = new mat_coef_t[num_materials];
    C2_[idx] = new mat_coef_t[num_materials];
    C3_[idx] = new mat_coef_t[num_materials];
    C4_[idx] = new mat_coef_t[num_materials];
    C5_[idx] = new mat_coef_t[num_materials];
    C6_[idx] = new mat_coef_t[num_materials];

    if (!Ca[idx]_ || !Cbx_[idx] || !Cbz_[idx] || !Cby_[idx] 
        || !Da_[idx] || !Dbx_[idx] || !Dbz_[idx] || !Dby_[idx]
        || !C1_[idx] || !C2_[idx] || !C3_[idx] || !C4_[idx] 
        || !C5_[idx] || !C6_[idx])
    {
      free_coefs();
      throw MemoryException();
    }
  }
}

void UPml::free_coefs()
{
  if (Ca_)
  {
    for (int idx = 0; idx < thicnkess_; idx++)
    {
      if (Ca_[idx])
        delete[] Ca_[idx];
      if (Cbx_[idx])
        delete[] Cbx_[idx];
      if (Cby_[idx])
        delete[] Cby_[idx];
      if (Cbz_[idx])
        delete[] Cbz_[idx];
      if (Da_[idx])
        delete[] Da_[idx];
      if (Dbx_[idx])
        delete[] Dbx_[idx];
      if (Dby_[idx])
        delete[] Dby_[idx];
      if (Dbz_[idx])
        delete[] Dbz_[idx];
      if (C1_[idx])
        delete[] C1_[idx];
      if (C2_[idx])
        delete[] C2_[idx];
      if (C3_[idx])
        delete[] C3_[idx];
      if (C4_[idx])
        delete[] C4_[idx];
      if (C5_[idx])
        delete[] C5_[idx];
      if (C6_[idx])
        delete[] C6_[idx];
    }
    delete[] Ca_;
    Ca_ = 0;
    delete[] Cbx_;
    delete[] Cby_;
    delete[] Cbz_;
    delete[] Da_;
    delete[] Dbx_;
    delete[] Dby_;
    delete[] Dbz_;
    delete[] C1_;
    delete[] C2_;
    delete[] C3_;
    delete[] C4_;
    delete[] C5_;
    delete[] C6_;
  }
}

void UPml::deinit(const Grid &grid, Face face)
{
  if (d_)
  {
    delete[] d_;
    d_ = 0;
  }

  if (h_)
  {
    delete[] h_;
    h_ = 0;
  }

  if (sigmas_)
  {
    delete[] sigmas_;
    sigmas_ = 0;
  }

  free_coefs();
}

void UPml::add_sd_bcs(SubdomainBc *sd, Face bcface, Face sdface)
{

}

void UPml::apply(Face face, Grid &grid, FieldType type)
{
  if (type == E)
  {
    switch (face)
    {
    case FRONT:
    case BACK:
      update_ex(grid, true);   
      update_ey(grid, false);
      update_ez(grid, false);
      break;

    case LEFT:
    case RIGHT:
      update_ex(grid, false);   
      update_ey(grid, true);
      update_ez(grid, false);
      break;

    case TOP:
    case BOTTOM:
      update_ex(grid, false);   
      update_ey(grid, false);
      update_ez(grid, true);
      break;
    }
  }
  else if (type == H)
  {
    switch (face)
    {
    case FRONT:
    case BACK:
      update_hx(grid, true);   
      update_hy(grid, false);
      update_hz(grid, false);
      break;

    case LEFT:
    case RIGHT:
      update_hx(grid, false);   
      update_hy(grid, true);
      update_hz(grid, false);
      break;

    case TOP:
    case BOTTOM:
      update_hx(grid, false);   
      update_hy(grid, false);
      update_hz(grid, true);
      break;
    }
  } 
  else
  {
    cout << "INCORRECT FIELD TYPE GIVEN TO UPDATE!" << endl;
  }
  
}

unsigned int UPml::sigma_index(Grid &grid, unsigned int x, 
                               unsigned int y, 
                               unsigned int z)
{
  unsigned int ret = 0;
  switch(our_face_)
  {
  case FRONT:
    ret = grid.get_ldx() - x;
    break;
  case BACK:
    ret = 
  case TOP:
  case BOTTOM:
  case LEFT:
  case RIGHT:
  }
}

void UPml::update_ex(Grid &grid, bool pml) 
{
  unsigned int grid_idx, pml_idx, mid, sig_idx = 0; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_ex_r_); 
  
  field_t d_temp = 0;

#ifdef USE_OPENMP
#pragam omp parallel private(mid, grid_idx, pml_idx, sig_idx, i, j, k, it, jt, kt, d_temp)
  {
#pragam omp for
#endif
    for(i = pml_r.xmin, it = grid_ex_r_.xmin; 
        it < grid_ex_r_.xmax; i++, it++, sig_idx++)
    {
      for(j = pml_r.ymin, jt = grid_ex_r_.ymin; 
          jt < grid_ex_r_.ymax; j++, jt++)
      {
        
        for(k = pml_r.zmin, kt = grid_ex_r_.zmin; 
            kt < grid_ex_r_.zmax; k++, kt++)
        {
          grid_idx = grid.pi(it, jt, kt);
          pml_idx = pi(i, j, k);
          
          mid = grid.material_[grid_idx];
          
          if (pml) 
          {
            d_temp = C1_[sig_idx][mid] * d_[pml_idx]
              + C2_[sig_idx][mid] 
              * (grid.hz_[grid.pi(it, jt, kt)] / grid.get_deltay() 
                 - grid.hz_[grid.pi(it, jt-1, kt)] / grid.get_deltaz());
            
            grid.ex_[grid_idx] = C3_[sig_idx][mid] * grid.ex_[grid_idx]
              + C4_[sig_idx][mid] * ( C5_[sig_idx][mid] * d_temp
                                     - C6_[sig_idx][mid] * d_[pml_idx]);
                                     
            d_[pml_idx] = d_temp;
                                     
          } else {
            grid.ex_[grid_idx] = Ca_[sig_idx][mid] * grid.ex_[grid_idx]
              + Cby_[sig_idx][mid] * (grid.hz_[grid.pi(it, jt, kt)] 
                                      - grid.hz_[grid.pi(it, jt-1, kt)])
              + Cbz_[sig_idx][mid] * (grid.hy_[grid_idx - 1] 
                                      - grid.hy_[grid_idx]);
          }
        }
      }
    }
#ifdef USE_OPENMP
  }
#endif
}

void UPml::update_ey(Grid &grid, bool pml) 
{
}

void UPml::update_ez(Grid &grid, bool pml) 
{
}

void UPml::update_hx(Grid &grid, bool pml)
{
}

void UPml::update_hy(Grid &grid, bool pml)
{
}

void UPml::update_hz(Grid &grid, bool pml)
{
}
