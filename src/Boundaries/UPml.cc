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

#include "UPml.hh"
#include "../Exceptions.hh"
#include "../Grid.hh"
#include "../Constants.hh"

#include "../config.h"

/** TEMP **/
/*#undef USE_OPENMP*/

#ifdef USE_OPENMP
#include <omp.h>
#endif

UPml::UPml()
  : common_(0), d_(0), h_(0)
{

}

UPml::~UPml()
{

}

void UPml::compute_regions(Face face, const Grid &grid)
{
  grid_r_ = find_face(face, grid);

  bc_r_.xmin = bc_r_.ymin = bc_r_.zmin = 0;

  bc_r_.xmax = grid_r_.xmax - grid_r_.xmin;
  bc_r_.ymax = grid_r_.ymax - grid_r_.ymin;
  bc_r_.zmax = grid_r_.zmax - grid_r_.zmin;

  grid_ex_r_ = grid_ey_r_ = grid_ez_r_ = grid_r_;
  grid_hx_r_ = grid_hy_r_ = grid_hz_r_ = grid_r_;

  grid_ex_r_.ymin++;
  grid_ex_r_.zmin++;
  
  grid_ey_r_.xmin++;
  grid_ey_r_.zmin++;
  
  grid_ez_r_.xmin++;
  grid_ez_r_.ymin++;
  
  grid_hx_r_.ymax--;
  grid_hx_r_.zmax--;
  
  grid_hy_r_.xmax--;
  grid_hy_r_.zmax--;
  
  grid_hz_r_.xmax--;
  grid_hz_r_.ymax--;

  // Corrections made by comparing to Jan's FDTD
  grid_ex_r_.xmax--;
  grid_ey_r_.ymax--;
  grid_ez_r_.zmax--;

  // Don't allow external face E field updates (electric walls)
  // Make sure that the PML computes all components at internal faces. 
  if (thickness_ > 0) 
  {
    switch (face) {
    case FRONT:
      grid_ey_r_.xmin--;
      grid_ez_r_.xmin--;

      grid_ey_r_.xmax--;
      grid_ez_r_.xmax--;

      grid_ex_r_.ymax--;
      grid_ez_r_.ymax--;

      grid_ex_r_.zmax--;
      grid_ey_r_.zmax--;
      break;

    case BACK:
      grid_ex_r_.ymax--;
      grid_ez_r_.ymax--;

      grid_ex_r_.zmax--;
      grid_ey_r_.zmax--;

      break;

    case TOP:
      grid_ex_r_.zmin--;
      grid_ey_r_.zmin--;

      grid_ex_r_.zmax--;
      grid_ey_r_.zmax--;

      grid_ey_r_.xmax--;
      grid_ez_r_.xmax--;

      grid_ez_r_.ymax--;
      grid_ex_r_.ymax--;
      break;

    case BOTTOM:
      grid_ey_r_.xmax--;
      grid_ez_r_.xmax--;

      grid_ez_r_.ymax--;
      grid_ex_r_.ymax--;
      break;

    case LEFT:
      grid_ey_r_.xmax--;
      grid_ez_r_.xmax--;

      grid_ex_r_.zmax--;
      grid_ey_r_.zmax--;     
      break;

    case RIGHT:
      grid_ez_r_.ymin--;
      grid_ex_r_.ymin--;

      grid_ez_r_.ymax--;
      grid_ex_r_.ymax--;

      grid_ey_r_.xmax--;
      grid_ez_r_.xmax--;

      grid_ex_r_.zmax--;
      grid_ey_r_.zmax--;
      break;
    }

    // Don't overlap compute corners and edges more than once
    // temporarially commented out...
//     switch (face)
//     {
//     case FRONT:
//     case BACK:
//       if (grid.get_boundary(LEFT).get_type() == UPML)
//       {
//         grid_ey_r_.ymin += grid.get_boundary(LEFT).get_thickness();
//         grid_hy_r_.ymin += grid.get_boundary(LEFT).get_thickness() + 1;
//       }

//       if (grid.get_boundary(RIGHT).get_type() == UPML)
//       {
//         grid_ey_r_.ymax -= grid.get_boundary(LEFT).get_thickness();
//         grid_hy_r_.ymax -= grid.get_boundary(LEFT).get_thickness() + 1;
//       }

//       if (grid.get_boundary(TOP).get_type() == UPML)
//       {
//         grid_ez_r_.zmax -= grid.get_boundary(TOP).get_thickness();
//         grid_hz_r_.zmax -= grid.get_boundary(TOP).get_thickness() + 1;
//       }

//       if (grid.get_boundary(BOTTOM).get_type() == UPML)
//       {
//         grid_ez_r_.zmin += grid.get_boundary(BOTTOM).get_thickness();
//         grid_hz_r_.zmin += grid.get_boundary(BOTTOM).get_thickness() + 1;
//       }
//       break;

//     case LEFT:
//     case RIGHT:
//       if (grid.get_boundary(FRONT).get_type() == UPML)
//       {
//         unsigned int thickness = grid.get_boundary(FRONT).get_thickness();
//         grid_ex_r_.xmax -= thickness;
//         grid_hx_r_.xmax -= thickness + 1;
//         grid_ez_r_.xmax -= thickness;
//         grid_hz_r_.xmax -= thickness + 1;
//       }

//       if (grid.get_boundary(BACK).get_type() == UPML)
//       {
//         unsigned int thickness = grid.get_boundary(BACK).get_thickness();
//         grid_ex_r_.xmin += thickness;
//         grid_hx_r_.xmin += thickness + 1;
//         grid_ez_r_.xmin += thickness;
//         grid_hz_r_.xmin += thickness + 1;
//       }

//       if (grid.get_boundary(TOP).get_type() == UPML)
//       {
//         unsigned int thickness = grid.get_boundary(TOP).get_thickness();
//         grid_ez_r_.zmax -= thickness;
//         grid_hz_r_.zmax -= thickness + 1;
//         grid_ex_r_.zmax -= thickness;
//         grid_hx_r_.zmax -= thickness + 1;
//       }

//       if (grid.get_boundary(BOTTOM).get_type() == UPML)
//       {
//         unsigned int thickness = grid.get_boundary(BOTTOM).get_thickness();
//         grid_ez_r_.zmin += thickness;
//         grid_hz_r_.zmin += thickness + 1;
//         grid_ex_r_.zmin += thickness;
//         grid_hx_r_.zmin += thickness + 1;
//       }
//       break;

//     case TOP:
//     case BOTTOM:
//       if (grid.get_boundary(LEFT).get_type() == UPML)
//       {
//         unsigned int thickness = grid.get_boundary(LEFT).get_thickness();
//         grid_ey_r_.ymin += thickness;
//         grid_hy_r_.ymin += thickness + 1;
//         grid_ex_r_.ymin += thickness;
//         grid_hx_r_.ymin += thickness + 1;
//       }

//       if (grid.get_boundary(RIGHT).get_type() == UPML)
//       {
//         unsigned int thickness = grid.get_boundary(LEFT).get_thickness();
//         grid_ey_r_.ymax -= thickness;
//         grid_hy_r_.ymax -= thickness + 1;
//         grid_ex_r_.ymax -= thickness;
//         grid_hx_r_.ymax -= thickness + 1;
//       }

//       if (grid.get_boundary(FRONT).get_type() == UPML)
//       {
//         unsigned int thickness = grid.get_boundary(FRONT).get_thickness();
//         grid_ex_r_.xmax -= thickness;
//         grid_hx_r_.xmax -= thickness + 1;
//         grid_ey_r_.xmax -= thickness;
//         grid_hy_r_.xmax -= thickness;
//       }

//       if (grid.get_boundary(BACK).get_type() == UPML)
//       {
//         unsigned int thickness = grid.get_boundary(BACK).get_thickness();
//         grid_ex_r_.xmin += thickness;
//         grid_hx_r_.xmin += thickness + 1;
//         grid_ey_r_.xmin += thickness;
//         grid_hy_r_.xmin += thickness;
//       }
//       break;
//     }
  }
}

void UPml::init(const Grid &grid, Face face)
{
  common_ = UPmlCommon::get_upml_common(const_cast<Grid &>(grid));

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

  cout << "UPML Update region for face " << face << ":"
       << "\n\tEx, x: " << grid_ex_r_.xmin << " -> " 
       << grid_ex_r_.xmax
       << ", y: " << grid_ex_r_.ymin << " -> " 
       << grid_ex_r_.ymax
       << ", z: " << grid_ex_r_.zmin << " -> " 
       << grid_ex_r_.zmax
       << "\n\tEy, x: " << grid_ey_r_.xmin << " -> " 
       << grid_ey_r_.xmax
       << ", y: " << grid_ey_r_.ymin << " -> " 
       << grid_ey_r_.ymax
       << ", z: " << grid_ey_r_.zmin << " -> " 
       << grid_ey_r_.zmax
       << "\n\tEz, x: " << grid_ez_r_.xmin << " -> " 
       << grid_ez_r_.xmax
       << ", y: " << grid_ez_r_.ymin << " -> " 
       << grid_ez_r_.ymax
       << ", z: " << grid_ez_r_.zmin << " -> " 
       << grid_ez_r_.zmax 
       << "\n\tHx, x: " << grid_hx_r_.xmin << " -> " 
       << grid_hx_r_.xmax
       << ", y: " << grid_hx_r_.ymin << " -> " 
       << grid_hx_r_.ymax
       << ", z: " << grid_hx_r_.zmin << " -> " 
       << grid_hx_r_.zmax
       << "\n\tHy, x: " << grid_hy_r_.xmin << " -> " 
       << grid_hy_r_.xmax
       << ", y: " << grid_hy_r_.ymin << " -> " 
       << grid_hy_r_.ymax
       << ", z: " << grid_hy_r_.zmin << " -> " 
       << grid_hy_r_.zmax
       << "\n\tHz, x: " << grid_hz_r_.xmin << " -> " 
       << grid_hz_r_.xmax
       << ", y: " << grid_hz_r_.ymin << " -> " 
       << grid_hz_r_.ymax
       << ", z: " << grid_hz_r_.zmin << " -> " 
       << grid_hz_r_.zmax << endl;

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

void UPml::update_ex(Grid &grid) 
{
  unsigned int grid_idx, pml_idx, mid, sig_idx = 0; 

  loop_idx_t i,j,k; 	/* indices in PML-layer */
  loop_idx_t it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_ex_r_); 

  field_t d_temp = 0;

  // Pointers
  const field_t *ex, *hz1, *hz2, *hy, *d;

  sig_idx = 0;

  // Why is this outside of the for loop? OpenMP thing? 
  i = pml_r.xmin;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, grid_idx, pml_idx, sig_idx, i, j, k, it, jt, kt, d_temp)
  {
#pragma omp for
#endif
    for(it = grid_ex_r_.xmin; 
        it < grid_ex_r_.xmax; i++, it++)
    {
      for(j = pml_r.ymin, jt = grid_ex_r_.ymin; 
          jt < grid_ex_r_.ymax; j++, jt++)
      {
        grid_idx = grid.pi(it, jt, kt);
        pml_idx = pi(i, j, k);
          
        ex = grid.get_pointer(grid_point(it,jt,kt), FC_EX);
        hz1 = grid.get_pointer(grid_point(it,jt,kt), FC_HZ);
        hz2 = grid.get_pointer(grid_point(it,jt-1,kt), FC_HZ);
        hy = grid.get_pointer(grid_point(it,jt,kt), FC_HY);

        d = &(d_[pml_idx]);

        for(k = pml_r.zmin, kt = grid_ex_r_.zmin; 
            kt < grid_ex_r_.zmax; k++, kt++)
        {
          mid = grid.material_[grid_idx];

          // Update equations go here!
          d_temp = *d * common_->Ay(j) 
            + common_->By(j) * ( (*hz1 - *hz2) - (*(hy - 1) - *hy) );

          *ex = *ex * common_->Az(k) 
            + common_->Bz(k) * common_->er(i,j,k) 
            * (d_temp * common_->Cx(i) - *d * common_->Dx(i));

          ex++;
          hz1++;
          hz2++;
          hy++;
        }
      }
    }
#ifdef USE_OPENMP
  }
#endif
}

void UPml::update_ey(Grid &grid) 
{
  unsigned int grid_idx, pml_idx, mid, sig_idx = 0; 

  int i,j,k; 	/* indices in PML-layer */
  int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_ey_r_); 
  
  field_t d_temp = 0;

  i = pml_r.xmin;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, grid_idx, pml_idx, sig_idx, i, j, k, it, jt, kt, d_temp)
  {
#pragma omp for
#endif
    for(it = grid_ey_r_.xmin; 
        it < grid_ey_r_.xmax; it++)
    {
      for(sig_idx = 0, j = pml_r.ymin, jt = grid_ey_r_.ymin; 
          jt < grid_ey_r_.ymax; j++, jt++, sig_idx++)
      {
        
        for(k = pml_r.zmin, kt = grid_ey_r_.zmin; 
            kt < grid_ey_r_.zmax; k++, kt++)
        {
          grid_idx = grid.pi(it, jt, kt);
          pml_idx = pi(i, j, k);
          
          mid = grid.material_[grid_idx];
          
          // Update equations go here!
        }
      }

      i++;
    }
#ifdef USE_OPENMP
  }
#endif
}

void UPml::update_ez(Grid &grid) 
{
  unsigned int grid_idx, pml_idx, mid, sig_idx = 0; 

  int i,j,k; 	/* indices in PML-layer */
  int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_ez_r_); 
  
  field_t d_temp = 0;

  i = pml_r.xmin;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, grid_idx, pml_idx, sig_idx, i, j, k, it, jt, kt, d_temp)
  {
#pragma omp for
#endif
    for(it = grid_ez_r_.xmin; 
        it < grid_ez_r_.xmax; it++)
    {
      for(j = pml_r.ymin, jt = grid_ez_r_.ymin; 
          jt < grid_ez_r_.ymax; j++, jt++)
      {
        
        for(sig_idx = 0, k = pml_r.zmin, kt = grid_ez_r_.zmin; 
            kt < grid_ez_r_.zmax; k++, kt++, sig_idx++)
        {
          grid_idx = grid.pi(it, jt, kt);
          pml_idx = pi(i, j, k);
          
          mid = grid.material_[grid_idx];
          
          // Update equations go here!
        }
      }
      i++;
    }
#ifdef USE_OPENMP
  }
#endif
}

void UPml::update_hx(Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid, sig_idx = 0; 

  int i,j,k; 	/* indices in PML-layer */
  int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_hx_r_); 
  
  field_t h_temp = 0;

  sig_idx = 0; i = pml_r.xmin;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, grid_idx, pml_idx, sig_idx, i, j, k, it, jt, kt, h_temp)
  {
#pragma omp for
#endif
    for(it = grid_hx_r_.xmin; 
        it < grid_hx_r_.xmax; it++)
    {
      for(j = pml_r.ymin, jt = grid_hx_r_.ymin; 
          jt < grid_hx_r_.ymax; j++, jt++)
      {
        
        for(k = pml_r.zmin, kt = grid_hx_r_.zmin; 
            kt < grid_hx_r_.zmax; k++, kt++)
        {
          grid_idx = grid.pi(it, jt, kt);
          pml_idx = pi(i, j, k);
          
          mid = grid.material_[grid_idx];
          
          // Update equations go here. 
        }
      }
      sig_idx++;
      i++;
    }
#ifdef USE_OPENMP
  }
#endif
}

void UPml::update_hy(Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid, sig_idx = 0; 

  int i,j,k; 	/* indices in PML-layer */
  int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_hy_r_); 
  
  field_t h_temp = 0;
  i = pml_r.xmin;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, grid_idx, pml_idx, sig_idx, i, j, k, it, jt, kt, h_temp)
  {
#pragma omp for
#endif
    for(it = grid_hy_r_.xmin; 
        it < grid_hy_r_.xmax; it++)
    {
      for(sig_idx = 0, j = pml_r.ymin, jt = grid_hy_r_.ymin; 
          jt < grid_hy_r_.ymax; j++, jt++, sig_idx++)
      {
        
        for(k = pml_r.zmin, kt = grid_hy_r_.zmin; 
            kt < grid_hy_r_.zmax; k++, kt++)
        {
          grid_idx = grid.pi(it, jt, kt);
          pml_idx = pi(i, j, k);
          
          mid = grid.material_[grid_idx];
          
          // Update equations go here. 
        }
      }
      i++;
    }
#ifdef USE_OPENMP
  }
#endif
}

void UPml::update_hz(Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid, sig_idx = 0; 

  int i,j,k; 	/* indices in PML-layer */
  int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_hz_r_); 
  
  field_t h_temp = 0;
  i = pml_r.xmin;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, grid_idx, pml_idx, sig_idx, i, j, k, it, jt, kt, h_temp)
  {
#pragma omp for
#endif
    for(it = grid_hz_r_.xmin; 
        it < grid_hz_r_.xmax; it++)
    {
      for(j = pml_r.ymin, jt = grid_hz_r_.ymin; 
          jt < grid_hz_r_.ymax; j++, jt++)
      {
        
        for(sig_idx = 0, k = pml_r.zmin, kt = grid_hz_r_.zmin; 
            kt < grid_hz_r_.zmax; k++, kt++, sig_idx++)
        {
          grid_idx = grid.pi(it, jt, kt);
          pml_idx = pi(i, j, k);
          
          mid = grid.material_[grid_idx];
          
          // Update equations go here. 

        }
      }
      i++;
    }
#ifdef USE_OPENMP
  }
#endif
}

void UPml::add_sd_bcs(SubdomainBc *sd, Face bcface, Face sdface)
{

}
