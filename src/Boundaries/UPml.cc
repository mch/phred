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
  : common_(0), dx_(0), dy_(0), dz_(0), 
    bx_(0), by_(0), bz_(0)
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
  
  // new gleefully throws exceptions on low memory.
  dx_ = new field_t[sz];
  dy_ = new field_t[sz];
  dz_ = new field_t[sz];

  bx_ = new field_t[sz];
  by_ = new field_t[sz];
  bz_ = new field_t[sz];

  memset(dx_, 0, sizeof(field_t) * sz);
  memset(dy_, 0, sizeof(field_t) * sz);
  memset(dz_, 0, sizeof(field_t) * sz);
  memset(bx_, 0, sizeof(field_t) * sz);
  memset(by_, 0, sizeof(field_t) * sz);
  memset(bz_, 0, sizeof(field_t) * sz);

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
  if (dx_)
  {
    delete[] dx_;
    delete[] dy_;
    delete[] dz_;
    dx_ = 0;
  }

  if (bx_)
  {
    delete[] bx_;
    delete[] by_;
    delete[] bz_;
    bx_ = 0;
  }
}

void UPml::apply(Face face, Grid &grid, FieldType type)
{
  if (type == E)
  {
    update_ex(grid);   
    update_ey(grid);
    update_ez(grid);
  }
  else if (type == H)
  {
    update_hx(grid);   
    update_hy(grid);
    update_hz(grid);
  } 
  else
  {
    cout << "INCORRECT FIELD TYPE GIVEN TO UPDATE!" << endl;
  }
  
}

void UPml::update_ex(Grid &grid) 
{
  unsigned int grid_idx, grid_idx2, pml_idx, mid; 

  loop_idx_t i,j,k; 	/* indices in PML-layer */
  loop_idx_t it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_ex_r_); 

  field_t d_temp = 0;

  // Pointers
  field_t *ex, *hz1, *hz2, *hy, *dx;

  // Inverted delta's
  field_t idy = 1 / grid.get_deltay();
  field_t idz = 1 / grid.get_deltaz();

  // This is outside of the for() statement so that OpenMP can
  // parallelize the loop.
  i = pml_r.xmin;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, grid_idx, pml_idx, i, j, k, it, jt, kt, d_temp)
  {
#pragma omp for
#endif
    for(it = grid_ex_r_.xmin; 
        it < grid_ex_r_.xmax; it++)
    {
      for(j = pml_r.ymin, jt = grid_ex_r_.ymin; 
          jt < grid_ex_r_.ymax; j++, jt++)
      {
        grid_idx = grid.pi(it, jt, grid_ex_r_.zmin);
        grid_idx2 = grid.pi(it, jt - 1, grid_ex_r_.zmin);
        pml_idx = pi(i, j, pml_r.zmin);
          
        ex = grid.get_ex_ptr(grid_idx);
        hz1 = grid.get_hz_ptr(grid_idx);
        hz2 = grid.get_hz_ptr(grid_idx2);
        hy = grid.get_hy_ptr(grid_idx);

        dx = &(dx_[pml_idx]);

        for(k = pml_r.zmin, kt = grid_ex_r_.zmin; 
            kt < grid_ex_r_.zmax; k++, kt++)
        {
          mid = grid.material_[grid_idx];

          // Update equations go here!
          d_temp = *dx * common_->Ay(jt) 
            + common_->By(jt) * ( idy*(*hz1 - *hz2) - idz*(*(hy - 1) - *hy) );

          *ex = *ex * common_->Az(kt) 
            + common_->Bz(kt) * common_->er(mid) 
            * (d_temp * common_->Cx(it) - *dx * common_->Dx(it));

          *dx = d_temp;

          dx++;
          ex++;
          hz1++;
          hz2++;
          hy++;
          grid_idx++;
        }
      }
      i++;
    }
#ifdef USE_OPENMP
  }
#endif
}

void UPml::update_ey(Grid &grid) 
{
  unsigned int grid_idx, pml_idx, mid; 

  int i,j,k; 	/* indices in PML-layer */
  int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_ey_r_); 
  
  field_t d_temp = 0;

  field_t *ey, *hx, *hz1, *hz2, *dy;

  // Inverted delta's
  field_t idx = 1 / grid.get_deltax();
  field_t idz = 1 / grid.get_deltaz();

  i = pml_r.xmin;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, grid_idx, pml_idx, i, j, k, it, jt, kt, d_temp)
  {
#pragma omp for
#endif
    for(it = grid_ey_r_.xmin; 
        it < grid_ey_r_.xmax; it++)
    {
      for(j = pml_r.ymin, jt = grid_ey_r_.ymin; 
          jt < grid_ey_r_.ymax; j++, jt++)
      {
        grid_idx = grid.pi(it, jt, grid_ey_r_.zmin);
        pml_idx = pi(i, j, pml_r.zmin);
          
        ey = grid.get_ey_ptr(grid_idx);
        hz1 = grid.get_hz_ptr(grid.pi(it-1, jt, grid_ey_r_.zmin));
        hz2 = grid.get_hz_ptr(grid_idx);
        hx = grid.get_hx_ptr(grid_idx);

        dy = &(dy_[pml_idx]);

        for(k = pml_r.zmin, kt = grid_ey_r_.zmin; 
            kt < grid_ey_r_.zmax; k++, kt++)
        {
          mid = grid.material_[grid_idx];

          if (jt == 9 && kt == 9)
          {
            cout << "Ey should have non zero hy, hz... it = " 
                 << it << ", hy[" << it << ", " << jt << ", " 
                 << kt << "] = " << grid.get_hy(it,jt,kt) << endl;
          }
                    
          // Update equations go here!
          d_temp = *dy * common_->Az(kt) 
            + common_->Bz(kt) * ( idz*(*hx - *(hx-1)) - idx*(*hz1 - *hz2));

          *ey = *ey * common_->Ax(it) 
            + common_->Bx(it) * common_->er(mid)
            * ( d_temp * common_->Cy(jt) - *dy * common_->Dy(jt) );

          *dy = d_temp;

          // Increment
          ey++;
          hz1++;
          hz2++;
          hx++;
          dy++;
          grid_idx++;
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
  unsigned int grid_idx, pml_idx, mid; 

  int i,j,k; 	/* indices in PML-layer */
  int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_ez_r_); 
  
  field_t d_temp = 0;
  
  field_t *ez, *hy1, *hy2, *hx1, *hx2, *dz;

  // Inverted delta's
  field_t idx = 1 / grid.get_deltax();
  field_t idy = 1 / grid.get_deltay();

  i = pml_r.xmin;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, grid_idx, pml_idx, i, j, k, it, jt, kt, d_temp)
  {
#pragma omp for
#endif
    for(it = grid_ez_r_.xmin; 
        it < grid_ez_r_.xmax; it++)
    {
      for(j = pml_r.ymin, jt = grid_ez_r_.ymin; 
          jt < grid_ez_r_.ymax; j++, jt++)
      {
        grid_idx = grid.pi(it, jt, grid_ez_r_.zmin);
        pml_idx = pi(i, j, pml_r.zmin);
          
        ez = grid.get_ez_ptr(grid_idx);

        hy1 = grid.get_hy_ptr(grid_idx);
        hy2 = grid.get_hy_ptr(grid.pi(it-1, jt, grid_ez_r_.zmin));
        hx1 = grid.get_hx_ptr(grid.pi(it, jt-1, grid_ez_r_.zmin));
        hx2 = grid.get_hx_ptr(grid_idx);

        dz = &(dz_[pml_idx]);

        for(k = pml_r.zmin, kt = grid_ez_r_.zmin; 
            kt < grid_ez_r_.zmax; k++, kt++)
        {
          mid = grid.material_[grid_idx];
          
          // Update equations go here!
          d_temp = *dz * common_->Ax(it)
            + common_->Bx(it) * ( idx*(*hy1 - *hy2) - idy*(*hx1 - *hx2) );

          *ez = *ez * common_->Ay(jt)
            + common_->By(jt) * common_->er(mid) 
            * (d_temp * common_->Cz(kt) - *dz * common_->Dz(kt));

          *dz = d_temp;

          ez++;
          hy1++;
          hy2++;
          hx1++;
          hx2++;
          dz++;
          grid_idx++;
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
  unsigned int grid_idx, grid_idx2, pml_idx, mid; 

  loop_idx_t i,j,k;   /* indices in PML-layer */
  loop_idx_t it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_hx_r_); 
  
  field_t b_temp = 0;
  field_t *hx, *ez1, *ez2, *ey, *bx;

  // Inverted delta's
  field_t idy = 1 / grid.get_deltay();
  field_t idz = 1 / grid.get_deltaz();

  i = pml_r.xmin;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, grid_idx, pml_idx, i, j, k, it, jt, kt, b_temp)
  {
#pragma omp for
#endif
    for(it = grid_hx_r_.xmin; 
        it < grid_hx_r_.xmax; it++)
    {
      for(j = pml_r.ymin, jt = grid_hx_r_.ymin; 
          jt < grid_hx_r_.ymax; j++, jt++)
      {
        grid_idx = grid.pi(it, jt, grid_hx_r_.zmin);
        grid_idx2 = grid.pi(it, jt+1, grid_hx_r_.zmin);

        pml_idx = pi(i, j, pml_r.zmin);
          
        hx = grid.get_hx_ptr(grid_idx);
        ey = grid.get_ey_ptr(grid_idx);
        ez1 = grid.get_ez_ptr(grid_idx);
        ez2 = grid.get_ez_ptr(grid_idx2);
        
        bx = &(bx_[pml_idx]);

        for(k = pml_r.zmin, kt = grid_hx_r_.zmin; 
            kt < grid_hx_r_.zmax; k++, kt++)
        {
          mid = grid.material_[grid_idx];
          
          // Update equations go here. 
          b_temp = *bx * common_->Ay(jt)
            + common_->By(jt) * ( idz*(*ey - *(ey + 1)) - idy*(*ez2 - *ez1) );

          *hx = *hx * common_->Az(kt)
            + common_->Bz(kt) * common_->ur(mid) 
            * ( b_temp * common_->Cx(it) - *bx * common_->Dx(it) );

          *bx = b_temp;

          hx++;
          ez1++;
          ez2++;
          ey++;
          bx++;
          grid_idx++;
        }
      }
      i++;
    }
#ifdef USE_OPENMP
  }
#endif
}

void UPml::update_hy(Grid &grid)
{
  unsigned int grid_idx, grid_idx2, pml_idx, mid; 

  int i,j,k; 	/* indices in PML-layer */
  int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_hy_r_); 
  
  field_t b_temp = 0;
  field_t *hy, *ex, *ez1, *ez2, *by;

  // Inverted delta's
  field_t idx = 1 / grid.get_deltax();
  field_t idz = 1 / grid.get_deltaz();

  i = pml_r.xmin;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, grid_idx, pml_idx, i, j, k, it, jt, kt, b_temp)
  {
#pragma omp for
#endif
    for(it = grid_hy_r_.xmin; 
        it < grid_hy_r_.xmax; it++)
    {
      for(j = pml_r.ymin, jt = grid_hy_r_.ymin; 
          jt < grid_hy_r_.ymax; j++, jt++)
      {
        grid_idx = grid.pi(it, jt, grid_hy_r_.zmin);
        grid_idx2 = grid.pi(it+1, jt, grid_hy_r_.zmin);

        pml_idx = pi(i, j, pml_r.zmin);
          
        hy = grid.get_hy_ptr(grid_idx);
        ex = grid.get_ex_ptr(grid_idx);
        ez1 = grid.get_ez_ptr(grid_idx2);
        ez2 = grid.get_ez_ptr(grid_idx);
        
        by = &(by_[pml_idx]);

        for(k = pml_r.zmin, kt = grid_hy_r_.zmin; 
            kt < grid_hy_r_.zmax; k++, kt++)
        {
          mid = grid.material_[grid_idx];

          // Update equations go here. 
          b_temp = *by * common_->Az(kt)
            + common_->Bz(kt) * ( idx*(*ez2 - *ez1) - idz*(*(ex + 1) - *ex) );

          *hy = *hy * common_->Ax(it) 
            + common_->Bx(it) * common_->ur(mid)
            * (b_temp * common_->Cy(jt) - *by * common_->Dy(jt));

          *by = b_temp;

          hy++;
          ex++;
          ez1++;
          ez2++;
          by++;
          grid_idx++;
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
  unsigned int grid_idx, grid_idx2, grid_idx3, pml_idx, mid; 

  int i,j,k; 	/* indices in PML-layer */
  int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_hz_r_); 
  
  field_t b_temp = 0;
  field_t *hz, *ey1, *ey2, *ex1, *ex2, *bz;

  // Inverted delta's
  field_t idx = 1 / grid.get_deltax();
  field_t idy = 1 / grid.get_deltay();

  i = pml_r.xmin;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, grid_idx, pml_idx, i, j, k, it, jt, kt, h_temp)
  {
#pragma omp for
#endif
    for(it = grid_hz_r_.xmin; 
        it < grid_hz_r_.xmax; it++)
    {
      for(j = pml_r.ymin, jt = grid_hz_r_.ymin; 
          jt < grid_hz_r_.ymax; j++, jt++)
      {
        grid_idx = grid.pi(it, jt, grid_hz_r_.zmin);
        grid_idx2 = grid.pi(it+1, jt, grid_hz_r_.zmin);
        grid_idx3 = grid.pi(it, jt+1, grid_hz_r_.zmin);

        pml_idx = pi(i, j, pml_r.zmin);
          
        hz = grid.get_hz_ptr(grid_idx);

        ex1 = grid.get_ex_ptr(grid_idx3);
        ex2 = grid.get_ex_ptr(grid_idx);
        ey1 = grid.get_ey_ptr(grid_idx);
        ey2 = grid.get_ey_ptr(grid_idx2);
        
        bz = &(bz_[pml_idx]);

        for(k = pml_r.zmin, kt = grid_hz_r_.zmin; 
            kt < grid_hz_r_.zmax; k++, kt++)
        {
          mid = grid.material_[grid_idx];
          
          // Update equations go here. 
          b_temp = *bz * common_->Ax(it)
            + common_->Bx(it) * ( idy*(*ex2 - *ex1) - idx*(*ey2 - *ey1) );

          *hz = *hz * common_->Ay(jt)
            + common_->By(jt) * common_->ur(mid)
            * ( b_temp * common_->Cz(kt) - *bz * common_->Dz(kt));

          *bz = b_temp;

          hz++;
          ex1++; ex2++;
          ey1++; ey2++;
          bz++;
          grid_idx++;
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
