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

UPml::UPml()
{

}

UPml::~UPml()
{

}

void UPml::init(const Grid &grid, Face face)
{
  compute_regions(face, grid);

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
}

void UPml::deinit(const Grid &grid, Face face)
{

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
      //upml_update_ex(grid);   
      update_ey(grid);
      update_ez(grid);
      break;

    case LEFT:
    case RIGHT:
      update_ex(grid);   
      //upml_update_ey(grid);
      update_ez(grid);
      break;

    case TOP:
    case BOTTOM:
      update_ex(grid);   
      update_ey(grid);
      //upml_update_ez(grid);
      break;
    }
  }
  else if (type == H)
  {
    switch (face)
    {
    case FRONT:
    case BACK:
      //upml_update_hx(grid);   
      update_hy(grid);
      update_hz(grid);
      break;

    case LEFT:
    case RIGHT:
      update_hx(grid);   
      //upml_update_hy(grid);
      update_hz(grid);
      break;

    case TOP:
    case BOTTOM:
      update_hx(grid);   
      update_hy(grid);
      //upml_update_hz(grid);
      break;
    }
  } 
  else
  {
    cout << "INCORRECT FIELD TYPE GIVEN TO UPDATE!" << endl;
  }
  
}

void UPml::update_ex(Grid &grid) 
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_ex_r_); 

#ifdef USE_OPENMP
#pragam omp parallel private(mid, grid_idx, pml_idx, i, j, k, it, jt, kt)
  {
#pragam omp for
#endif
    for(i = pml_r.xmin, it = grid_ex_r_.xmin; it < grid_ex_r_.xmax; i++, it++)
      for(j = pml_r.ymin, jt = grid_ex_r_.ymin; jt < grid_ex_r_.ymax; j++, jt++)
        for(k = pml_r.zmin, kt = grid_ex_r_.zmin; kt < grid_ex_r_.zmax; k++, kt++)
        {
          grid_idx = grid.pi(it, jt, kt);
          pml_idx = pi(i, j, k);
          
          mid = grid.material_[grid_idx];
          
          d_[pml_idx] = d_[pml_idx] + grid.get_dt() 
            * ((grid.hy_[grid.pi(it, jt, kt-1)] 
                - grid.hy_[grid_idx]) / grid.get_dz() 
               - (grid.hz_[grid_idx] 
                  - grid.hz_[grid.pi(it, jt-1, kt)]) / grid.get_dy());

          grid.ex_[grid_idx] = grid.ex_[grid_idx] + 1/(EPS_0 * 
        }
#ifdef USE_OPENMP
  }
#endif
}

void UPml::update_ey(Grid &grid) 
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
          
          *ey = grid.Ca_[mid] * *ey
            + grid.Cbz_[mid] * (*hx - *(hx-1))
            + grid.Cbx_[mid] * (*hz1 - *hz2);
          
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

void UPml::update_ez(Grid &grid) 
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
          
          *ez = grid.Ca_[mid] * *ez
            + grid.Cbx_[mid] * (*hy1 - *hy2)
            + grid.Cby_[mid] * (*hx1 - *hx2);

          ez++;
          hy1++; hy2++; hx1++; hx2++;
          idx++;
        }
      }
    }
  }
}

void UPml::update_hx(Grid &grid)
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
          
          *hx = grid.Da_[mid] * *hx
            + grid.Dby_[mid] * (*ez1 - *ez2)
            + grid.Dbz_[mid] * (*(ey+1) - *ey);
          
          hx++; idx++;
          ez1++; ez2++; ey++;
        }
      }
    }
  }

}

void UPml::update_hy(Grid &grid)
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
          
          *hy = grid.Da_[mid] * *hy
            + grid.Dbz_[mid] * (*ex - *(ex + 1))
            + grid.Dbx_[mid] * (*ez1 - *ez2);        
          
          hy++; idx++;
          ex++; ez1++; ez2++;
        }
      }
    }
  }
}

void UPml::update_hz(Grid &grid)
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
          
          *hz1 = grid.Da_[mid] * *hz1
            + grid.Dbx_[mid] * (*ey1 - *ey2)
            + grid.Dby_[mid] * (*ex1 - *ex2);
          
          hz1++; idx++;
          ey1++; ey2++;
          ex1++; ex2++;
        }
      }
    }
  }
}
