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
#include "../Types.hh"

/** TEMP **/
/*#undef USE_OPENMP*/

#ifdef USE_OPENMP
#include <omp.h>
#endif

UPml::UPml()
  : common_(0), dx_(0), dy_(0), dz_(0), 
    bx_(0), by_(0), bz_(0), 
    aux1_x_(0), aux1_y_(0), aux1_z_(0),
    aux2_x_(0), aux2_y_(0), aux2_z_(0),
    aux3_x_(0), aux3_y_(0), aux3_z_(0),
    sigma_max_(0.0), 
    poly_order_(4), eps_opt_(1.0), sigma_ratio_(1.0)
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

    case TOP:
      grid_ex_r_.zmin--;
      grid_ey_r_.zmin--;

      grid_ex_r_.zmax--;
      grid_ey_r_.zmax--;

      break;

    case RIGHT:
      grid_ez_r_.ymin--;
      grid_ex_r_.ymin--;

      grid_ez_r_.ymax--;
      grid_ex_r_.ymax--;

      break;
    }

    // Don't overlap compute corners and edges more than once
    if (face == LEFT || face == RIGHT || face == TOP || face == BOTTOM)
    {
      // Don't overlap the edges with the front and the back
      if (grid.get_boundary(BACK).get_type() == UPML)
      {
        unsigned int thickness = grid.get_boundary(BACK).get_thickness();
        grid_ex_r_.xmin += thickness;
        grid_ey_r_.xmin += thickness;
        grid_ez_r_.xmin += thickness;

        grid_hx_r_.xmin += thickness + 1;
        grid_hy_r_.xmin += thickness;
        grid_hz_r_.xmin += thickness;
      }

      // Could this block of code be causing the BOTTOM and FRONT
      // problem when divided along the Z axis?
      if (grid.get_boundary(FRONT).get_type() == UPML)
      {
        unsigned int thickness = grid.get_boundary(FRONT).get_thickness();
        
        grid_ex_r_.xmax -= thickness;
        grid_ey_r_.xmax -= thickness;
        grid_ez_r_.xmax -= thickness;

        grid_hx_r_.xmax -= thickness + 1;
        grid_hy_r_.xmax -= thickness;
        grid_hz_r_.xmax -= thickness;
      }      
    }

    if (face == TOP || face == BOTTOM)
    {
      if (grid.get_boundary(LEFT).get_type() == UPML)
      {
        unsigned int thickness = grid.get_boundary(LEFT).get_thickness();
        
        grid_ex_r_.ymin += thickness;
        grid_ey_r_.ymin += thickness;
        grid_ez_r_.ymin += thickness;

        grid_hx_r_.ymin += thickness;
        grid_hy_r_.ymin += thickness + 1;
        grid_hz_r_.ymin += thickness;
      }

      if (grid.get_boundary(RIGHT).get_type() == UPML)
      {
        unsigned int thickness = grid.get_boundary(RIGHT).get_thickness();
        
        grid_ex_r_.ymax -= thickness;
        grid_ey_r_.ymax -= thickness;
        grid_ez_r_.ymax -= thickness;

        grid_hx_r_.ymax -= thickness;
        grid_hy_r_.ymax -= thickness + 1;
        grid_hz_r_.ymax -= thickness;
      }
    }

    // Clean up, make E walls on outer edges
    switch (face) {

    case FRONT:
    case BACK:
      if (grid.get_boundary(RIGHT).get_type() == UPML)
      {
        grid_ex_r_.ymax--;
        grid_ez_r_.ymax--;
      }

      if (grid.get_boundary(TOP).get_type() == UPML)
      {
        grid_ex_r_.zmax--;
        grid_ey_r_.zmax--;
      }
      break;

    case TOP:
    case BOTTOM:
      if (grid.get_boundary(FRONT).get_type() == UPML)
      {
        grid_ey_r_.xmax--;
        grid_ez_r_.xmax--;
      }

      if (grid.get_boundary(RIGHT).get_type() == UPML)
      {
        grid_ez_r_.ymax--;
        grid_ex_r_.ymax--;
      }
      break;
      
    case LEFT:
    case RIGHT:
      if (grid.get_boundary(FRONT).get_type() == UPML)
      {
        grid_ey_r_.xmax--;
        grid_ez_r_.xmax--;
      }

      if (grid.get_boundary(TOP).get_type() == UPML)
      {
        grid_ex_r_.zmax--;
        grid_ey_r_.zmax--;
      }
      break;
    }

  }
}

void UPml::init(const Grid &grid, Face face)
{
  common_ = UPmlCommon::get_upml_common(const_cast<Grid &>(grid));

  common_->set_sigma_max(sigma_max_);
  common_->set_poly_order(poly_order_);
  common_->set_eps_opt(eps_opt_);
  common_->set_sigma_ratio(sigma_ratio_);

  common_->init_coeffs();

  compute_regions(face, grid);

  // A bit of a waste since some can be eliminated because there is no
  // overlap.
  unsigned int sz = (grid_r_.xmax - grid_r_.xmin + 1) 
    * (grid_r_.ymax - grid_r_.ymin + 1) 
    * (grid_r_.zmax - grid_r_.zmin + 1);
  
  // new gleefully throws exceptions on low memory.
  dx_ = new field_t[sz];
  dy_ = new field_t[sz];
  dz_ = new field_t[sz];

  bx_ = new field_t[sz];
  by_ = new field_t[sz];
  bz_ = new field_t[sz];

  cout << "UPml::init(). sz = " << sz << " field_t's. \ndx_ = " 
       << dx_ << " -> " << dx_ + sz << "\ndy_ = " << dy_ 
       << " -> " << dy_ + sz << "\ndz_ = " << dz_ << " -> " 
       << dz_ + sz << "\nbx_ = " 
       << bx_ << " -> " << bx_ + sz << "\nby_ = " << by_ << " -> " 
       << by_ + sz << "\nbz_ = " << bz_ << " -> " << bz_ + sz << endl;

  memset(dx_, 0, sizeof(field_t) * sz);
  memset(dy_, 0, sizeof(field_t) * sz);
  memset(dz_, 0, sizeof(field_t) * sz);
  memset(bx_, 0, sizeof(field_t) * sz);
  memset(by_, 0, sizeof(field_t) * sz);
  memset(bz_, 0, sizeof(field_t) * sz);

  shared_ptr<MaterialLib> mlib = grid.get_material_lib();
  map<string, Material>::const_iterator iter 
    = mlib->get_material_iter_begin();
  map<string, Material>::const_iterator iter_e 
    = mlib->get_material_iter_end();
  
  bool need_lossy = false;
  bool need_drude = false;
  bool need_debye = false;

  while (iter != iter_e) 
  {
    if ((*iter).second.type() == LOSSY)
      need_lossy = true;

    if ((*iter).second.type() == DRUDE)
      need_drude = true;

    if ((*iter).second.type() == DEBYE)
      need_debye = true;

    ++iter;
  }

  if (need_lossy || need_drude || need_debye)
  {
    aux1_x_ = new field_t[sz];
    aux1_y_ = new field_t[sz];
    aux1_z_ = new field_t[sz];
    
    memset(aux1_x_, 0, sizeof(field_t) * sz);
    memset(aux1_y_, 0, sizeof(field_t) * sz);
    memset(aux1_z_, 0, sizeof(field_t) * sz);
  }

  if (need_drude)
  {
    aux2_x_ = new field_t[sz];
    aux2_y_ = new field_t[sz];
    aux2_z_ = new field_t[sz];
    
    memset(aux2_x_, 0, sizeof(field_t) * sz);
    memset(aux2_y_, 0, sizeof(field_t) * sz);
    memset(aux2_z_, 0, sizeof(field_t) * sz);

    aux3_x_ = new field_t[sz];
    aux3_y_ = new field_t[sz];
    aux3_z_ = new field_t[sz];
    
    memset(aux3_x_, 0, sizeof(field_t) * sz);
    memset(aux3_y_, 0, sizeof(field_t) * sz);
    memset(aux3_z_, 0, sizeof(field_t) * sz);

  }

#ifdef DEBUG
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


  cout << endl << "br_r_: " << bc_r_;
  cout << "grid_r_: " << grid_r_ << endl;
#endif

  MPI_Type_contiguous(bc_r_.zmax, GRID_MPI_TYPE, &z_vector_);
  MPI_Type_commit(&z_vector_);

  cout << "UPml::init(): z vector is " << bc_r_.zmax << " long. " << endl;

  MPI_Type_vector(bc_r_.ymax, 1, bc_r_.zmax, GRID_MPI_TYPE, &y_vector_);
  MPI_Type_commit(&y_vector_);

  cout << "UPml::init(): y vector is " << bc_r_.ymax 
       << " long, stride is " << bc_r_.zmax << endl;
  
  MPI_Type_vector(bc_r_.xmax, 1, bc_r_.ymax * bc_r_.zmax, 
                  GRID_MPI_TYPE, &x_vector_);
  MPI_Type_commit(&x_vector_);

  cout << "UPml::init(): x vector is " << bc_r_.xmax 
       << " long, stride is " << bc_r_.ymax * bc_r_.zmax << endl;

  MPI_Type_contiguous(bc_r_.zmax * bc_r_.ymax, GRID_MPI_TYPE, &yz_plane_);
  MPI_Type_commit(&yz_plane_);

  cout << "UPml::init(): yz plane is a contiguous chunk of memory, " 
       << bc_r_.zmax * bc_r_.ymax << " total cells, " << bc_r_.zmax 
       << " by " << bc_r_.ymax << endl;

  MPI_Type_vector(bc_r_.xmax, bc_r_.zmax, bc_r_.ymax * bc_r_.zmax, 
                  GRID_MPI_TYPE, &xz_plane_);
  MPI_Type_commit(&xz_plane_);

  cout << "UPml::init(): xz plane is " << bc_r_.xmax << " by " << bc_r_.zmax 
       << ", stride is " << bc_r_.ymax * bc_r_.zmax << endl;


  // SUSPECT!
  MPI_Type_hvector(bc_r_.xmax, 1, sizeof(field_t) * bc_r_.zmax * bc_r_.ymax, 
                   y_vector_, &xy_plane_);
//   MPI_Type_vector(bc_r_.xmax, 1, bc_r_.zmax, 
//                   y_vector_, &xy_plane_);
  MPI_Type_commit(&xy_plane_);

  cout << "UPml::init(): xy plane is " << bc_r_.xmax << " by " << bc_r_.ymax 
       << ", stride is " << bc_r_.ymax * bc_r_.zmax << endl;
  
}

void UPml::deinit(const Grid &grid, Face face)
{
  if (common_)
    common_->deinit();

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

  if (aux1_x_)
  {
    delete[] aux1_x_;
    delete[] aux1_y_;
    delete[] aux1_z_;

    aux1_x_ = 0;
    aux1_y_ = 0;
    aux1_z_ = 0;
  }

  if (aux2_x_)
  {
    delete[] aux2_x_;
    delete[] aux2_y_;
    delete[] aux2_z_;

    aux2_x_ = 0;
    aux2_y_ = 0;
    aux2_z_ = 0;
  }

  if (aux3_x_)
  {
    delete[] aux3_x_;
    delete[] aux3_y_;
    delete[] aux3_z_;

    aux3_x_ = 0;
    aux3_y_ = 0;
    aux3_z_ = 0;
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
  loop_idx_t jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_ex_r_); 

  field_t d_temp = 0;

  // Pointers
  field_t *ex, *hz1, *hz2, *hy, *dx;

  // Inverted delta's
  field_t idy = 1 / grid.get_deltay();
  field_t idz = 1 / grid.get_deltaz();
  field_t idt = 1 / grid.get_deltat();

  // This is outside of the for() statement so that OpenMP can
  // parallelize the loop.
  int offset = grid_ex_r_.xmin - pml_r.xmin;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, grid_idx, pml_idx, i, j, k, jt, kt, \
                             d_temp, grid_idx2, ex, hz1, hz2, hy, dx)
  {
#pragma omp for
#endif
    for(loop_idx_t it = grid_ex_r_.xmin; 
        it < grid_ex_r_.xmax; it++)
    {
      i = it - offset;

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
          if (common_->mtype(mid) == PERF_COND) {
            *ex = *ex; 
          } 
          else if (common_->mtype(mid) == LOSSY) 
          {
            field_t q_temp;
            q_temp = aux1_x_[pml_idx] * common_->lossy_A(mid)
              + common_->lossy_B(mid) 
              * ( idy*(*hz1 - *hz2) - idz*(*hy - *(hy - 1)) );

            d_temp = common_->Ay(jt) * *dx
              + idt * common_->By(jt) * (q_temp - aux1_x_[pml_idx]);

            aux1_x_[pml_idx] = q_temp;

            *ex = *ex * common_->Az(kt) 
              + common_->Bz(kt) * common_->er(mid) 
              * (d_temp * common_->Cx(it) - *dx * common_->Dx(it));
            
            *dx = d_temp;
          } 
          else if (common_->mtype(mid) == DRUDE) 
          {
            // Unmagnatized plasma update

            // ADE fields:
            // p_temp: current value of P, P^n
            // aux1: previous value of P, P^(n-1)

            // d_temp: current value of D, D^n
            // *d: previous value of D, D^(n-1)
            // aux2: penultimate value of D, D^(n-2)

            // e_temp: current value of E, E^n
            // *e: previous value of E, E^(n-1)
            // aux3: penultimate value of E, E^(n-2)

            // Z transform fields:
            // aux2: Previous value of S, S^(n-1)
            // aux3: Penultimate value of S, S^(n-2)

            field_t p_temp; 
            p_temp = aux1_x_[pml_idx] * common_->Ay(jt) 
              + common_->By(jt) 
              * ( idy*(*hz1 - *hz2) - idz*(*hy - *(hy - 1)) );

            d_temp = *dx * common_->Az(kt) 
            + common_->Bz(kt)
            * (p_temp * common_->Cx(it) - aux1_x_[pml_idx] * common_->Dx(it));

            // E field update using ADE method
#ifdef ADE_DRUDE
            field_t e_temp;
            e_temp = common_->drude_c1(mid) * *ex 
              + common_->drude_c2(mid) * aux3_x_[pml_idx]
              + common_->drude_c3(mid) * d_temp
              + common_->drude_c4(mid) * *dx
              + common_->drude_c5(mid) * aux2_x_[pml_idx]
              ;

            // Advance storage locations. 
            aux2_x_[pml_idx] = *dx;
            *dx = d_temp;

            aux3_x_[pml_idx] = *ex;
            *ex = e_temp;
#else
            // E field update using Z transform method... Have to use
            // the same as the Grid does, which is Z transform. Using
            // ADE in UPML results in instability due to the small
            // error.
            float s_temp;
            *ex = d_temp - aux2_x_[pml_idx];
            
            s_temp = 
              (1 + common_->get_vcdt(mid)) * aux2_x_[pml_idx]
              - (common_->get_vcdt(mid) * aux3_x_[pml_idx])
              + (common_->get_omegasq(mid) * (1 - common_->get_vcdt(mid)))
              * *ex;

            // Advance storage locations. 
            aux3_x_[pml_idx] = aux2_x_[pml_idx];
            aux2_x_[pml_idx] = s_temp;
#endif

            // Advance storage locations. 
            aux1_x_[pml_idx] = p_temp;
          }
          else if (common_->mtype(mid) == DEBYE) 
          {
            field_t p_temp; 
            p_temp = aux1_x_[pml_idx] * common_->Ay(jt) 
              + common_->By(jt) 
              * ( idy*(*hz1 - *hz2) - idz*(*hy - *(hy - 1)) );

            d_temp = *dx * common_->Az(kt) 
            + common_->Bz(kt)
            * (p_temp * common_->Cx(it) - aux1_x_[pml_idx] * common_->Dx(it));

            *ex = *ex * common_->debyeA(mid) 
              + d_temp * common_->debyeB(mid)
              - *dx * common_->debyeC(mid);

            *dx = d_temp;
            aux1_x_[pml_idx] = p_temp;
          }
          else 
          { // DIELECTRIC
            d_temp = *dx * common_->Ay(jt) 
              + common_->By(jt) 
              * ( idy*(*hz1 - *hz2) - idz*(*hy - *(hy - 1)) );

            *ex = *ex * common_->Az(kt) 
              + common_->Bz(kt) * common_->er(mid) 
              * (d_temp * common_->Cx(it) - *dx * common_->Dx(it));

            *dx = d_temp;
          }

          dx++;
          ex++;
          hz1++;
          hz2++;
          hy++;
          grid_idx++;
          pml_idx++;
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

  loop_idx_t i,j,k; 	/* indices in PML-layer */
  loop_idx_t jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_ey_r_); 
  
  field_t d_temp = 0;

  field_t *ey, *hx, *hz1, *hz2, *dy;

  // Inverted delta's
  field_t idx = 1 / grid.get_deltax();
  field_t idz = 1 / grid.get_deltaz();
  field_t idt = 1 / grid.get_deltat();

  loop_idx_t offset = grid_ey_r_.xmin - pml_r.xmin;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, grid_idx, pml_idx, i, j, k, jt, kt, \
                             d_temp, ey, hx, hz1, hz2, dy)
  {
#pragma omp for
#endif
    for(loop_idx_t it = grid_ey_r_.xmin; 
        it < grid_ey_r_.xmax; it++)
    {
      i = it - offset;

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

          // Update equations go here!
          if (common_->mtype(mid) == PERF_COND) {
            *ey = *ey;
          } 
          else if (common_->mtype(mid) == LOSSY) 
          {
            field_t q_temp;
            q_temp = aux1_y_[pml_idx] * common_->lossy_A(mid)
              + common_->lossy_B(mid) 
              * ( idz*(*hx - *(hx-1)) - idx*(*hz2 - *hz1));

            d_temp = common_->Az(kt) * *dy
              + idt * common_->Bz(kt) * (q_temp - aux1_y_[pml_idx]);

            aux1_y_[pml_idx] = q_temp;

            *ey = *ey * common_->Ax(it) 
              + common_->Bx(it) * common_->er(mid)
              * ( d_temp * common_->Cy(jt) - *dy * common_->Dy(jt) );
            
            *dy = d_temp;
          } 
          else if (common_->mtype(mid) == DRUDE) 
          {
            // Unmagnatized plasma update

            // ADE fields:
            // p_temp: current value of P, P^n
            // aux1: previous value of P, P^(n-1)

            // d_temp: current value of D, D^n
            // *d: previous value of D, D^(n-1)
            // aux2: penultimate value of D, D^(n-2)

            // e_temp: current value of E, E^n
            // *e: previous value of E, E^(n-1)
            // aux3: penultimate value of E, E^(n-2)

            // Z transform fields:
            // aux2: Previous value of S, S^(n-1)
            // aux3: Penultimate value of S, S^(n-2)

            field_t p_temp; 
            p_temp = aux1_y_[pml_idx] * common_->Az(jt) 
              + common_->Bz(jt) 
              * ( idz*(*hx - *(hx-1)) - idx*(*hz2 - *hz1));

            d_temp = *dy * common_->Ax(kt) 
            + common_->Bx(kt)
            * (p_temp * common_->Cy(it) - aux1_y_[pml_idx] * common_->Dy(it));

#ifdef ADE_DRUDE
            field_t e_temp;
            e_temp = common_->drude_c1(mid) * *ey 
              + common_->drude_c2(mid) * aux3_y_[pml_idx]
              + common_->drude_c3(mid) * d_temp
              + common_->drude_c4(mid) * *dy
              + common_->drude_c5(mid) * aux2_y_[pml_idx]
              ;

            aux2_y_[pml_idx] = *dy;
            *dy = d_temp;

            aux3_y_[pml_idx] = *ey;
            *ey = e_temp;
#else
            // E field update using Z transform method... Have to use
            // the same as the Grid does, which is Z transform. Using
            // ADE in UPML results in instability due to the small
            // error.
            float s_temp;
            *ey = d_temp - aux2_y_[pml_idx];
            
            s_temp = 
              (1 + common_->get_vcdt(mid)) * aux2_y_[pml_idx]
              - (common_->get_vcdt(mid) * aux3_y_[pml_idx])
              + (common_->get_omegasq(mid) * (1 - common_->get_vcdt(mid)))
              * *ey;

            // Advance storage locations. 
            aux3_y_[pml_idx] = aux2_y_[pml_idx];
            aux2_y_[pml_idx] = s_temp;
#endif

            // Advance storage locations. 
            aux1_y_[pml_idx] = p_temp;

          }
          else if (common_->mtype(mid) == DEBYE) 
          {
            field_t p_temp; 
            p_temp = aux1_y_[pml_idx] * common_->Az(kt) 
              + common_->Bz(kt) 
              * ( idz*(*hx - *(hx-1)) - idx*(*hz2 - *hz1));

            d_temp = *dy * common_->Ax(it) 
            + common_->Bx(it)
            * (p_temp * common_->Cy(jt) - aux1_y_[pml_idx] * common_->Dy(jt));

            *ey = *ey * common_->debyeA(mid) 
              + d_temp * common_->debyeB(mid)
              - *dy * common_->debyeC(mid);

            *dy = d_temp;
            aux1_y_[pml_idx] = p_temp;
          }
          else 
          { // DIELECTRIC
            d_temp = *dy * common_->Az(kt) 
              + common_->Bz(kt) 
              * ( idz*(*hx - *(hx-1)) - idx*(*hz2 - *hz1));

            *ey = *ey * common_->Ax(it) 
              + common_->Bx(it) * common_->er(mid)
              * ( d_temp * common_->Cy(jt) - *dy * common_->Dy(jt) );

            *dy = d_temp;
          }

          // Increment
          ey++;
          hz1++;
          hz2++;
          hx++;
          dy++;
          grid_idx++;
          pml_idx++;
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

  loop_idx_t i,j,k; 	/* indices in PML-layer */
  loop_idx_t jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_ez_r_); 
  
  field_t d_temp = 0;
  
  field_t *ez, *hy1, *hy2, *hx1, *hx2, *dz;

  // Inverted delta's
  field_t idx = 1 / grid.get_deltax();
  field_t idy = 1 / grid.get_deltay();
  field_t idt = 1 / grid.get_deltat();

  loop_idx_t offset = grid_ez_r_.xmin - pml_r.xmin;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, grid_idx, pml_idx, i, j, k, jt, kt, \
                             d_temp, ez, hy1, hy2, hx1, hx2, dz)
  {
#pragma omp for
#endif
    for(loop_idx_t it = grid_ez_r_.xmin; 
        it < grid_ez_r_.xmax; it++)
    {
      i = it - offset;

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
          if (common_->mtype(mid) == PERF_COND) {
            *ez = *ez;
          } 
          else if (common_->mtype(mid) == LOSSY) 
          {
            field_t q_temp;
            q_temp = aux1_z_[pml_idx] * common_->lossy_A(mid)
              + common_->lossy_B(mid) 
              * ( idx*(*hy1 - *hy2) - idy*(*hx2 - *hx1) );

            d_temp = common_->Ax(it) * *dz
              + idt * common_->Bx(it) * (q_temp - aux1_z_[pml_idx]);

            aux1_z_[pml_idx] = q_temp;

            *ez = *ez * common_->Ay(jt)
              + common_->By(jt) * common_->er(mid) 
              * (d_temp * common_->Cz(kt) - *dz * common_->Dz(kt));

            *dz = d_temp;
          } 
          else if (common_->mtype(mid) == DRUDE) 
          {
            // Unmagnatized plasma update

            // ADE fields:
            // p_temp: current value of P, P^n
            // aux1: previous value of P, P^(n-1)

            // d_temp: current value of D, D^n
            // *d: previous value of D, D^(n-1)
            // aux2: penultimate value of D, D^(n-2)

            // e_temp: current value of E, E^n
            // *e: previous value of E, E^(n-1)
            // aux3: penultimate value of E, E^(n-2)

            // Z transform fields:
            // aux2: Previous value of S, S^(n-1)
            // aux3: Penultimate value of S, S^(n-2)

            field_t p_temp; 
            p_temp = aux1_z_[pml_idx] * common_->Ax(jt) 
              + common_->Bx(jt) 
              * ( idx*(*hy1 - *hy2) - idy*(*hx2 - *hx1) );

            d_temp = *dz * common_->Ay(kt) 
            + common_->By(kt)
            * (p_temp * common_->Cz(it) - aux1_z_[pml_idx] * common_->Dz(it));

#ifdef ADE_DRUDE
            e_temp = common_->drude_c1(mid) * *ez 
              + common_->drude_c2(mid) * aux3_z_[pml_idx]
              + common_->drude_c3(mid) * d_temp
              + common_->drude_c4(mid) * *dz
              + common_->drude_c5(mid) * aux2_z_[pml_idx]
              ;

            aux2_z_[pml_idx] = *dz;
            *dz = d_temp;

            aux3_z_[pml_idx] = *ez;
            *ez = e_temp;
#else
            // E field update using Z transform method... Have to use
            // the same as the Grid does, which is Z transform. Using
            // ADE in UPML results in instability due to the small
            // error.
            float s_temp;
            *ez = d_temp - aux2_z_[pml_idx];
            
            s_temp = 
              (1 + common_->get_vcdt(mid)) * aux2_z_[pml_idx]
              - (common_->get_vcdt(mid) * aux3_z_[pml_idx])
              + (common_->get_omegasq(mid) * (1 - common_->get_vcdt(mid)))
              * *ez;

            // Advance storage locations. 
            aux3_z_[pml_idx] = aux2_z_[pml_idx];
            aux2_z_[pml_idx] = s_temp;
#endif
            // Advance storage locations. 
            aux1_z_[pml_idx] = p_temp;

          }
          else if (common_->mtype(mid) == DEBYE) 
          {
            field_t p_temp; 
            p_temp = aux1_z_[pml_idx] * common_->Ax(it) 
              + common_->Bx(it) 
              * ( idx*(*hy1 - *hy2) - idy*(*hx2 - *hx1) );

            d_temp = *dz * common_->Ay(jt) 
            + common_->By(jt)
            * (p_temp * common_->Cz(kt) - aux1_y_[pml_idx] * common_->Dz(kt));

            *ez = *ez * common_->debyeA(mid) 
              + d_temp * common_->debyeB(mid)
              - *dz * common_->debyeC(mid);

            *dz = d_temp;
            aux1_z_[pml_idx] = p_temp;
          }
          else 
          { // DIELECTRIC
            d_temp = *dz * common_->Ax(it)
              + common_->Bx(it) 
              * ( idx*(*hy1 - *hy2) - idy*(*hx2 - *hx1) );

            *ez = *ez * common_->Ay(jt)
              + common_->By(jt) * common_->er(mid) 
              * (d_temp * common_->Cz(kt) - *dz * common_->Dz(kt));

            *dz = d_temp;
          }

          ez++;
          hy1++;
          hy2++;
          hx1++;
          hx2++;
          dz++;
          grid_idx++;
          pml_idx++;
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
  loop_idx_t jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_hx_r_); 
  
  field_t b_temp = 0;
  field_t *hx, *ez1, *ez2, *ey, *bx;

  // Inverted delta's
  field_t idy = 1 / grid.get_deltay();
  field_t idz = 1 / grid.get_deltaz();

  loop_idx_t offset = grid_hx_r_.xmin - pml_r.xmin;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, grid_idx, pml_idx, i, j, k, jt, kt, \
                             b_temp, hx, ez1, ez2, ey, bx, grid_idx2)
  {
#pragma omp for
#endif
    for(loop_idx_t it = grid_hx_r_.xmin; 
        it < grid_hx_r_.xmax; it++)
    {
      i = it - offset;

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
          if (common_->mtype(mid) == PERF_COND) {
            *hx = *hx; 
          } 
          else 
          {
            b_temp = *bx * common_->Ay(jt)
              + common_->By(jt) 
              * ( idz*(*(ey + 1) - *ey) - idy*(*ez2 - *ez1) );

            *hx = *hx * common_->Az(kt)
              + common_->Bz(kt) * common_->ur(mid) 
              * ( b_temp * common_->Cx(it) - *bx * common_->Dx(it) );
            
            *bx = b_temp;
          }

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

  loop_idx_t i,j,k; 	/* indices in PML-layer */
  loop_idx_t jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_hy_r_); 
  
  field_t b_temp = 0;
  field_t *hy, *ex, *ez1, *ez2, *by;

  // Inverted delta's
  field_t idx = 1 / grid.get_deltax();
  field_t idz = 1 / grid.get_deltaz();

  loop_idx_t offset = grid_hy_r_.xmin - pml_r.xmin;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, grid_idx, pml_idx, i, j, k, jt, kt, \
                             b_temp, grid_idx2, hy, ex, ez1, ez2, by)
  {
#pragma omp for
#endif
    for(loop_idx_t it = grid_hy_r_.xmin; 
        it < grid_hy_r_.xmax; it++)
    {
      i = it - offset;

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
          if (common_->mtype(mid) == PERF_COND) {
            *hy = *hy; 
          } 
          else 
          {
            b_temp = *by * common_->Az(kt)
              + common_->Bz(kt) 
              * ( idx*(*ez1 - *ez2) - idz*(*(ex + 1) - *ex) );

            *hy = *hy * common_->Ax(it) 
              + common_->Bx(it) * common_->ur(mid)
              * (b_temp * common_->Cy(jt) - *by * common_->Dy(jt));

            *by = b_temp;
          }

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

  loop_idx_t i,j,k; 	/* indices in PML-layer */
  loop_idx_t jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_hz_r_); 
  
  field_t b_temp = 0;
  field_t *hz, *ey1, *ey2, *ex1, *ex2, *bz;

  // Inverted delta's
  field_t idx = 1 / grid.get_deltax();
  field_t idy = 1 / grid.get_deltay();

  loop_idx_t offset = grid_hz_r_.xmin - pml_r.xmin;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, grid_idx, pml_idx, i, j, k, jt, kt, \
                             b_temp, grid_idx2, grid_idx3, hz, ey1, ey2, \
                             ex1, ex2, bz)
  {
#pragma omp for
#endif
    for(loop_idx_t it = grid_hz_r_.xmin; 
        it < grid_hz_r_.xmax; it++)
    {
      i = it - offset;

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
          if (common_->mtype(mid) == PERF_COND) {
            *hz = *hz; 
          } 
          else 
          {
            b_temp = *bz * common_->Ax(it)
              + common_->Bx(it) 
              * ( idy*(*ex1 - *ex2) - idx*(*ey2 - *ey1) );

            *hz = *hz * common_->Ay(jt)
              + common_->By(jt) * common_->ur(mid)
              * ( b_temp * common_->Cz(kt) - *bz * common_->Dz(kt));

            *bz = b_temp;
          }

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
  RxTxData ret;
  unsigned int s_idx = 0, r_idx = 0;

  switch (sdface) {
  case BACK:
    s_idx = pi(1, 0, 0);
    r_idx = pi(0, 0, 0);
    ret.set_datatype(yz_plane_);
    break;

  case FRONT:
    s_idx = pi(bc_r_.xmax - 2, 0, 0);
    r_idx = pi(bc_r_.xmax - 1, 0, 0);
    ret.set_datatype(yz_plane_);
    break;

  case LEFT:
    s_idx = pi(0, 1, 0);
    r_idx = pi(0, 0, 0);
    ret.set_datatype(xz_plane_);
    break;

  case RIGHT:
    s_idx = pi(0, bc_r_.ymax - 2, 0);
    r_idx = pi(0, bc_r_.ymax - 1, 0);
    ret.set_datatype(xz_plane_);
    break;

  case BOTTOM:
    s_idx = pi(0, 0, 1);
    r_idx = pi(0, 0, 0);
    ret.set_datatype(xy_plane_);
    break;

  case TOP:
    s_idx = pi(0, 0, bc_r_.zmax - 2);
    r_idx = pi(0, 0, bc_r_.zmax - 1);
    ret.set_datatype(xy_plane_);
    break;
  }

  ret.set_field_type(E);

  ret.set_rx_ptr(&(dx_[r_idx]));
  ret.set_tx_ptr(&(dx_[s_idx]));
  sd->add_tx_rx_data(ret);

//   cout << "UPml::add_sd_bcs, bcface = " << face_string(bcface) 
//        << ", sdface = " << face_string(sdface) 
//        << ", dx, rx_ptr = " << &(dx_[r_idx]) << ", sx_ptr = "
//        << &(dx_[s_idx]) << endl;

  ret.set_rx_ptr(&(dy_[r_idx]));
  ret.set_tx_ptr(&(dy_[s_idx]));
  sd->add_tx_rx_data(ret);

//   cout << "UPml::add_sd_bcs, dy, rx_ptr = " << &(dy_[r_idx]) << ", sx_ptr = "
//        << &(dy_[s_idx]) << endl;

  ret.set_rx_ptr(&(dz_[r_idx]));
  ret.set_tx_ptr(&(dz_[s_idx]));
  sd->add_tx_rx_data(ret);

//   cout << "UPml::add_sd_bcs, dz, rx_ptr = " << &(dz_[r_idx]) << ", sx_ptr = "
//        << &(dz_[s_idx]) << endl;

  ret.set_field_type(H);

  ret.set_rx_ptr(&(bx_[r_idx]));
  ret.set_tx_ptr(&(bx_[s_idx]));
  sd->add_tx_rx_data(ret);

  ret.set_rx_ptr(&(by_[r_idx]));
  ret.set_tx_ptr(&(by_[s_idx]));
  sd->add_tx_rx_data(ret);

  ret.set_rx_ptr(&(bz_[r_idx]));
  ret.set_tx_ptr(&(bz_[s_idx]));
  sd->add_tx_rx_data(ret);

  // All aux variables are used in the computation of the E field. 
  ret.set_field_type(E);

  if (aux1_x_)
  {
    ret.set_rx_ptr(&(aux1_x_[r_idx]));
    ret.set_tx_ptr(&(aux1_x_[s_idx]));
    sd->add_tx_rx_data(ret);

    ret.set_rx_ptr(&(aux1_y_[r_idx]));
    ret.set_tx_ptr(&(aux1_y_[s_idx]));
    sd->add_tx_rx_data(ret);

    ret.set_rx_ptr(&(aux1_z_[r_idx]));
    ret.set_tx_ptr(&(aux1_z_[s_idx]));
    sd->add_tx_rx_data(ret);
  }

  if (aux2_x_)
  {
    ret.set_rx_ptr(&(aux2_x_[r_idx]));
    ret.set_tx_ptr(&(aux2_x_[s_idx]));
    sd->add_tx_rx_data(ret);

    ret.set_rx_ptr(&(aux2_y_[r_idx]));
    ret.set_tx_ptr(&(aux2_y_[s_idx]));
    sd->add_tx_rx_data(ret);

    ret.set_rx_ptr(&(aux2_z_[r_idx]));
    ret.set_tx_ptr(&(aux2_z_[s_idx]));
    sd->add_tx_rx_data(ret);
  }

  if (aux3_x_)
  {
    ret.set_rx_ptr(&(aux3_x_[r_idx]));
    ret.set_tx_ptr(&(aux3_x_[s_idx]));
    sd->add_tx_rx_data(ret);

    ret.set_rx_ptr(&(aux3_y_[r_idx]));
    ret.set_tx_ptr(&(aux3_y_[s_idx]));
    sd->add_tx_rx_data(ret);

    ret.set_rx_ptr(&(aux3_z_[r_idx]));
    ret.set_tx_ptr(&(aux3_z_[s_idx]));
    sd->add_tx_rx_data(ret);
  }

}
