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

  // Allocate memory
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

  for (int idx = 0; idx < thickness_; idx++)
    sigmas_[idx] = sigma_max * pow(static_cast<float>(idx), 
                                   static_cast<float>(poly_order_)) 
      / pow(static_cast<float>(thickness_),
            static_cast<float>(poly_order_));

  // Build the shadow library
  
  
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

void UPml::update_ex(Grid &grid, bool pml) 
{
  unsigned int grid_idx, pml_idx, mid, sig_idx; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_ex_r_); 
  
  field_t d_temp = 0;

#ifdef USE_OPENMP
#pragam omp parallel private(mid, grid_idx, pml_idx, sig_idx, i, j, k, it, jt, kt)
  {
#pragam omp for
#endif
    for(i = pml_r.xmin, it = grid_ex_r_.xmin; it < grid_ex_r_.xmax; i++, it++)
    {
      for(j = pml_r.ymin, jt = grid_ex_r_.ymin; jt < grid_ex_r_.ymax; j++, jt++)
      {

        ex = &(ex_[idx]);
        hz1 = &(hz_[idx]);
        hz2 = &(hz_[idx2]);
        hy = &(hy_[idx]);

        for(k = pml_r.zmin, kt = grid_ex_r_.zmin; kt < grid_ex_r_.zmax; k++, kt++)
        {
          grid_idx = grid.pi(it, jt, kt);
          pml_idx = pi(i, j, k);
          
          mid = grid.material_[grid_idx];
          
          if (pml) 
          {
            d_temp = d_[pml_idx] + grid.get_delta() 
              * ((grid.hy_[grid.pi(it, jt, kt-1)] 
                  - grid.hy_[grid_idx]) / grid.get_delta() 
                 - (grid.hz_[grid_idx] 
                    - grid.hz_[grid.pi(it, jt-1, kt)]) / grid.get_delta());

            grid.ex_[grid_idx] = grid.ex_[grid_idx] + 1/(EPS_0) 
              * (d_temp * (1 + (sigmas_[sig_idx] * grid.get_delta())/(2 * EPS_0))
                 - d_[pml_idx] * (1 - (sigmas_[sig_idx] * grid.get_delta())/(2 * EPS_0)));

            d_[pml_idx] = d_temp;

          } else {
            *ex = Ca_[mid] * *ex
              + Cby_[mid] * (*hz1 - *hz2)
              + Cbz_[mid] * (*(hy - 1) - *hy);
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
  unsigned int mid, idx;
  int i, j, k;
  field_t *ey, *hx, *hz1, *hz2;

  // Inner part
#ifdef USE_OPENMP
#pragma omp parallel private(mid, i, j, k, idx, ey, hx, hz1, hz2)
#endif
  {
#ifdef USE_OPENMP
#pragma omp for
#endif
    for (i = grid_ey_r_.xmin; i < grid_ey_r_.xmax; i++) {
      for (j = grid_ey_r_.ymin; j < grid_ey_r_.ymax; j++) {
        
        idx = pi(i, j, grid_ey_r_.zmin);
        hz1 = &(grid.hz_[pi(i-1, j, grid_ey_r_.zmin)]);

        hz2 = &(grid.hz_[idx]);
        ey = &(grid.ey_[idx]);
        hx = &(grid.hx_[idx]);

        for (k = grid_ey_r_.zmin; k < grid_ey_r_.zmax; k++) {
          mid = grid.material_[idx];
          
          if (pml)
          {

          } else {
            *ey = Ca_[mid] * *ey
              + Cbz_[mid] * (*hx - *(hx-1))
              + Cbx_[mid] * (*hz1 - *hz2);
          }
            
          ey++;
          hx++;
          hz1++;
          hz2++;
          idx++;
        }
      }
    }
  }
}

void UPml::update_ez(Grid &grid, bool pml) 
{
  unsigned int mid, idx;
  int i, j, k;
  field_t *ez, *hy1, *hy2, *hx1, *hx2;
  
  // Inner part
#ifdef USE_OPENMP
#pragma omp parallel private(mid, i, j, k, idx, ez, hy1, hy2, hx1, hx2)
#endif
  {
#ifdef USE_OPENMP
#pragma omp for
#endif
    for (i = grid_ez_r_.xmin; i < grid_ez_r_.xmax; i++) {
      for (j = grid_ez_r_.ymin; j < grid_ez_r_.ymax; j++) {

        idx = pi(i, j, grid_ez_r_.zmin);
        hy2 = &(grid.hy_[pi(i-1, j, grid_ez_r_.zmin)]);
        hx1 = &(grid.hx_[pi(i, j-1, grid_ez_r_.zmin)]);

        hx2 = &(grid.hx_[idx]);
        ez = &(grid.ez_[idx]);
        hy1 = &(grid.hy_[idx]);

        for (k = grid_ez_r_.zmin; k < grid_ez_r_.zmax; k++) {
          mid = grid.material_[idx];
          
          if (pml)
          {

          } else {
            *ez = Ca_[mid] * *ez
              + Cbx_[mid] * (*hy1 - *hy2)
              + Cby_[mid] * (*hx1 - *hx2);
          }

          ez++;
          hy1++; hy2++; hx1++; hx2++;
          idx++;
        }
      }
    }
  }
}

void UPml::update_hx(Grid &grid, bool pml)
{
  unsigned int mid, idx;
  int i, j, k;
  field_t *hx, *ez1, *ez2, *ey;

#ifdef USE_OPENMP
#pragma omp parallel  private(mid, i, j, k, idx, hx, ez1, ez2, ey)
#endif
  {
#ifdef USE_OPENMP
#pragma omp for
#endif
    for (i = grid_hx_r_.xmin; i < grid_hx_r_.xmax; i++) {
      for (j = grid_hx_r_.ymin; j < grid_hx_r_.ymax; j++) {
        
        idx = pi(i, j, grid_hx_r_.zmin);
        ez2 = &(grid.ez_[pi(i, j+1, grid_hx_r_.zmin)]);

        ey = &(grid.ey_[idx]);
        hx = &(grid.hx_[idx]);
        ez1 = &(grid.ez_[idx]);

        for (k = grid_hx_r_.zmin; k < grid_hx_r_.zmax; k++) {
          mid = grid.material_[idx];
          
          if (pml)
          {

          } else {
            *hx = Da_[mid] * *hx
              + Dby_[mid] * (*ez1 - *ez2)
              + Dbz_[mid] * (*(ey+1) - *ey);
          }

          hx++; idx++;
          ez1++; ez2++; ey++;
        }
      }
    }
  }

}

void UPml::update_hy(Grid &grid, bool pml)
{
  unsigned int mid, idx;
  int i, j, k;
  field_t *hy, *ex, *ez1, *ez2;


#ifdef USE_OPENMP
#pragma omp parallel private(mid, i, j, k, idx, hy, ex, ez1, ez2)
#endif
  {
#ifdef USE_OPENMP
#pragma omp for
#endif
    for (i = grid_hy_r_.xmin; i < grid_hy_r_.xmax; i++) {
      for (j = grid_hy_r_.ymin; j < grid_hy_r_.ymax; j++) {

        idx = pi(i, j, grid_hy_r_.zmin);
        ez1 = &(grid.ez_[pi(i+1, j, grid_hy_r_.zmin)]);

        hy = &(grid.hy_[idx]);
        ex = &(grid.ex_[idx]);
        ez2 = &(grid.ez_[idx]);

        for (k = grid_hy_r_.zmin; k < grid_hy_r_.zmax; k++) {
          mid = grid.material_[idx];
          
          if (pml)
          {

          } else {
            *hy = Da_[mid] * *hy
              + Dbz_[mid] * (*ex - *(ex + 1))
              + Dbx_[mid] * (*ez1 - *ez2);        
          }

          hy++; idx++;
          ex++; ez1++; ez2++;
        }
      }
    }
  }
}

void UPml::update_hz(Grid &grid, bool pml)
{
  unsigned int mid, idx;
  int i, j, k;
  field_t *hz1, *ey1, *ey2, *ex1, *ex2;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, i, j, k, idx, hz1, ey1, ey2, ex1, ex2)
#endif
  {
#ifdef USE_OPENMP
#pragma omp for
#endif
    for (i = grid_hz_r_.xmin; i < grid_hz_r_.xmax; i++) {
      for (j = grid_hz_r_.ymin; j < grid_hz_r_.ymax; j++) {
        
        idx = pi(i, j, grid_hz_r_.zmin);
        ey2 = &(grid.ey_[pi(i+1, j, grid_hz_r_.zmin)]);
        ex1 = &(grid.ex_[pi(i, j+1, grid_hz_r_.zmin)]);

        ex2 = &(grid.ex_[idx]);
        hz1 = &(grid.hz_[idx]);
        ey1 = &(grid.ey_[idx]);

        for (k = grid_hz_r_.zmin; k < grid_hz_r_.zmax; k++) {
          mid = grid.material_[idx];
          
          if (pml)
          {

          } else {
            *hz1 = Da_[mid] * *hz1
              + Dbx_[mid] * (*ey1 - *ey2)
              + Dby_[mid] * (*ex1 - *ex2);
          }

          hz1++; idx++;
          ey1++; ey2++;
          ex1++; ex2++;
        }
      }
    }
  }
}
