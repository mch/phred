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

#include "Pml.hh"
#include "Grid.hh"
#include "Exceptions.hh"
#include "Constants.hh"

#include <math.h>

using namespace std;

Pml::Pml()
  : variation_(VP), g_(0.0), nrml_refl_(1.0),
    ratio_m_(0.0), exponent_n_(0.0), delta_bndy_(0.0), 
    geometric_delta_(0.0), geometric_profile_(0),
    exy_(0), exz_(0), eyx_(0), eyz_(0), ezx_(0), ezy_(0),
    hxy_(0), hxz_(0), hyx_(0), hyz_(0), hzx_(0), hzy_(0), alloced_(false)
{}

Pml::Pml(PmlVariation_t variation, float g, float nrml_refl)
  : variation_(variation), g_(g), nrml_refl_(nrml_refl),
    ratio_m_(0.0), exponent_n_(0.0), delta_bndy_(0.0), 
    geometric_delta_(0.0), geometric_profile_(0),
    exy_(0), exz_(0), eyx_(0), eyz_(0), ezx_(0), ezy_(0),
    hxy_(0), hxz_(0), hyx_(0), hyz_(0), hzx_(0), hzy_(0), alloced_(false)
{}

Pml::Pml(PmlVariation_t variation, float nrml_refl)
  : variation_(variation), g_(0.0), nrml_refl_(nrml_refl),
    ratio_m_(0.0), exponent_n_(0.0), delta_bndy_(0.0), 
    geometric_delta_(0.0), geometric_profile_(0),
    exy_(0), exz_(0), eyx_(0), eyz_(0), ezx_(0), ezy_(0),
    hxy_(0), hxz_(0), hyx_(0), hyz_(0), hzx_(0), hzy_(0), alloced_(false)
{}

Pml::~Pml()
{
  free_pml_fields();
}

Pml::Pml(const Pml &rhs)
{
  *this = rhs;
}

const Pml &Pml::operator=(const Pml &rhs)
{
  variation_ = rhs.variation_;
  g_ = rhs.g_;
  nrml_refl_ = rhs.nrml_refl_;
  ratio_m_ = rhs.ratio_m_;
  exponent_n_ = rhs.exponent_n_;
  delta_bndy_ = rhs.delta_bndy_;
  geometric_delta_ = rhs.geometric_delta_;
  geometric_profile_ = rhs.geometric_profile_;
  exy_ = exz_ = eyx_ = eyz_ = ezx_ = ezy_ = 0;
  hxy_ = hxz_ = hyx_ = hyz_ = hzx_ = hzy_ = 0;
  alloced_ = false;

  thickness_ = rhs.thickness_;

  return *this;
}

void Pml::alloc_pml_fields(Face face, Grid &grid)
{
  if (alloced_)
    return; 

  if (thickness_ == 0)
    throw exception();
  
  region_t r = find_face(face, grid);

  pml_r_.xmin = pml_r_.ymin = pml_r_.zmin = 0;

  pml_r_.xmax = r.xmax - r.xmin;
  pml_r_.ymax = r.ymax - r.ymin;
  pml_r_.zmax = r.zmax - r.zmin;

  grid_ex_r_ = grid_ey_r_ = grid_ez_r_ = r;
  grid_hx_r_ = grid_hy_r_ = grid_hz_r_ = r;

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
      //grid_hz_r_.xmax++;
      //grid_hy_r_.xmax++;

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
      //grid_hy_r_.zmax++;
      //grid_hx_r_.zmax++;

      grid_ey_r_.xmax--;
      grid_ez_r_.xmax--;

      grid_ez_r_.ymax--;
      grid_ex_r_.ymax--;
      break;

    case LEFT:
      //grid_hz_r_.ymax++;
      //grid_hx_r_.ymax++;

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
  }

//   cout << "PML Update region for face " << face << ":"
//        << "\n\tEx, x: " << grid_ex_r_.xmin << " -> " 
//        << grid_ex_r_.xmax
//        << ", y: " << grid_ex_r_.ymin << " -> " 
//        << grid_ex_r_.ymax
//        << ", z: " << grid_ex_r_.zmin << " -> " 
//        << grid_ex_r_.zmax
//        << "\n\tEy, x: " << grid_ey_r_.xmin << " -> " 
//        << grid_ey_r_.xmax
//        << ", y: " << grid_ey_r_.ymin << " -> " 
//        << grid_ey_r_.ymax
//        << ", z: " << grid_ey_r_.zmin << " -> " 
//        << grid_ey_r_.zmax
//        << "\n\tEz, x: " << grid_ez_r_.xmin << " -> " 
//        << grid_ez_r_.xmax
//        << ", y: " << grid_ez_r_.ymin << " -> " 
//        << grid_ez_r_.ymax
//        << ", z: " << grid_ez_r_.zmin << " -> " 
//        << grid_ez_r_.zmax 
//        << "\n\tHx, x: " << grid_hx_r_.xmin << " -> " 
//        << grid_hx_r_.xmax
//        << ", y: " << grid_hx_r_.ymin << " -> " 
//        << grid_hx_r_.ymax
//        << ", z: " << grid_hx_r_.zmin << " -> " 
//        << grid_hx_r_.zmax
//        << "\n\tHy, x: " << grid_hy_r_.xmin << " -> " 
//        << grid_hy_r_.xmax
//        << ", y: " << grid_hy_r_.ymin << " -> " 
//        << grid_hy_r_.ymax
//        << ", z: " << grid_hy_r_.zmin << " -> " 
//        << grid_hy_r_.zmax
//        << "\n\tHz, x: " << grid_hz_r_.xmin << " -> " 
//        << grid_hz_r_.xmax
//        << ", y: " << grid_hz_r_.ymin << " -> " 
//        << grid_hz_r_.ymax
//        << ", z: " << grid_hz_r_.zmin << " -> " 
//        << grid_hz_r_.zmax << endl;

  unsigned int sz = (r.xmax - r.xmin) * (r.ymax - r.ymin) 
    * (r.zmax - r.zmin);
  
  exy_ = new field_t[sz];
  exz_ = new field_t[sz];

  eyx_ = new field_t[sz];
  eyz_ = new field_t[sz];

  ezx_ = new field_t[sz];
  ezy_ = new field_t[sz];

  hxy_ = new field_t[sz];
  hxz_ = new field_t[sz];

  hyx_ = new field_t[sz];
  hyz_ = new field_t[sz];

  hzx_ = new field_t[sz];
  hzy_ = new field_t[sz];
  
  if (exy_ && exz_ && eyx_ && eyz_ && ezx_ && ezy_
      && exy_ && exz_ && eyx_ && eyz_ && ezx_ && ezy_)
  {
    alloced_ = true;
  } else {
    free_pml_fields();
    throw MemoryException(); // Insufficent memory
  }

  memset(exy_, 0, sizeof(field_t) * sz);  
  memset(exz_, 0, sizeof(field_t) * sz);  
  memset(eyx_, 0, sizeof(field_t) * sz);  
  memset(eyz_, 0, sizeof(field_t) * sz);  
  memset(ezx_, 0, sizeof(field_t) * sz);  
  memset(ezy_, 0, sizeof(field_t) * sz);  

  memset(hxy_, 0, sizeof(field_t) * sz);  
  memset(hxz_, 0, sizeof(field_t) * sz);  
  memset(hyx_, 0, sizeof(field_t) * sz);  
  memset(hyz_, 0, sizeof(field_t) * sz);  
  memset(hzx_, 0, sizeof(field_t) * sz);  
  memset(hzy_, 0, sizeof(field_t) * sz);  

}

void Pml::set_thickness(unsigned int thickness)
{
  thickness_ = thickness;
}

void Pml::setup(Face face, Grid &grid)
{
  delta_t d_space;

  switch(face)
  {
  case FRONT:
  case BACK:
    d_space = grid.get_deltax();
    break;

  case TOP:
  case BOTTOM:
    d_space = grid.get_deltaz();
    break;

  case LEFT:
  case RIGHT:
    d_space = grid.get_deltay();
  }

  ratio_m_ = (thickness_ == 0 || nrml_refl_ > 9.0) ? 0.0 :
    (-log(nrml_refl_ * 0.01)) * C / (2.0 * thickness_ * d_space);

  delta_bndy_ = (thickness_ == 0) ? d_space : thickness_ * d_space;

  geometric_delta_ = d_space;

  geometric_profile_ = 0;
  switch (variation_) 
  {
  case VC:
    exponent_n_ = 0.0;
    break;

  case VL:
    exponent_n_ = 1.0;
    ratio_m_ *= 2;
    break;

  case VP:
    exponent_n_ = 2.0;
    ratio_m_ *= 3;
    break;

  case VG:
    ratio_m_ *= log(g_) / (pow(g_, static_cast<float>(thickness_)) - 1.0);
    geometric_profile_ = 1;
    break;
  }

  alloc_pml_fields(face, grid);

  MPI_Datatype y_vec;

  MPI_Type_vector(pml_r_.ymax, 1, pml_r_.zmax, GRID_MPI_TYPE, &y_vec);

  MPI_Type_contiguous(pml_r_.zmax * pml_r_.ymax, GRID_MPI_TYPE, &yz_plane_);
  MPI_Type_commit(&yz_plane_);

  MPI_Type_vector(pml_r_.xmax, pml_r_.zmax, pml_r_.ymax * pml_r_.zmax, 
                  GRID_MPI_TYPE, &xz_plane_);
  MPI_Type_commit(&xz_plane_);

  MPI_Type_hvector(pml_r_.xmax, 1, sizeof(field_t) * pml_r_.ymax 
                   * pml_r_.zmax, y_vec, &xy_plane_);
  MPI_Type_commit(&xy_plane_);
}

void Pml::free_pml_fields()
{
  if (exy_) 
    delete[] exy_;

  if (exz_) 
    delete[] exz_;

  if (eyx_) 
    delete[] eyx_;

  if (eyz_) 
    delete[] eyz_;

  if (ezx_) 
    delete[] ezx_;

  if (ezy_) 
    delete[] ezy_;

  if (hxy_) 
    delete[] hxy_;

  if (hxz_) 
    delete[] hxz_;

  if (hyx_) 
    delete[] hyx_;

  if (hyz_) 
    delete[] hyz_;

  if (hzx_) 
    delete[] hzx_;

  if (hzy_) 
    delete[] hzy_;

  exy_ = exz_ = eyx_ = eyz_ = ezx_ = ezy_ = hxy_ = hxz_ = 0;
  hyx_ = hyz_ = hzx_ = hzy_ = 0;
}

void Pml::apply(Face face, Grid &grid, FieldType type)
{
  if (!alloced_)
    throw exception(); // PML must be set up before applying it. 

  region_t grid_r = find_face(face, grid);

  if (type == E)
  {
    pml_update_ex(grid);   
    pml_update_ey(grid);
    pml_update_ez(grid);

  }
  else if (type == H)
  {
    pml_update_hx(grid);
    pml_update_hy(grid);
    pml_update_hz(grid);
  } 
  else
  {
    cout << "INCORRECT FIELD TYPE GIVEN TO UPDATE!" << endl;
  }
  
}

void Pml::pml_update_ex(Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  region_t pml_r; // Region in the PML to update

  pml_r.xmin = pml_r.ymin = pml_r.zmin = 0;

  pml_r.xmax = grid_ex_r_.xmax - grid_ex_r_.xmin;
  pml_r.ymax = grid_ex_r_.ymax - grid_ex_r_.ymin;
  pml_r.zmax = grid_ex_r_.zmax - grid_ex_r_.zmin;

  PmlCommon &com = grid.get_pml_common();

  for(i = pml_r.xmin, it = grid_ex_r_.xmin; it < grid_ex_r_.xmax; i++, it++)
    for(j = pml_r.ymin, jt = grid_ex_r_.ymin; jt < grid_ex_r_.ymax; j++, jt++)
      for(k = pml_r.zmin, kt = grid_ex_r_.zmin; kt < grid_ex_r_.zmax; k++, kt++)
      {
        grid_idx = grid.pi(it, jt, kt);
        pml_idx = pi(i, j, k);

        mid = grid.material_[grid_idx];

        exz_[pml_idx] = 
          com.get_e_z_coef1(kt) * grid.Ca_[mid] * exz_[pml_idx] 
          + com.get_e_z_coef2(kt) * grid.Cbz_[mid] 
            * (grid.hy_[grid.pi(it, jt, kt-1)] 
               - grid.hy_[grid.pi(it, jt, kt)]);
        
        exy_[pml_idx] = 
          com.get_e_y_coef1(jt) * grid.Ca_[mid] * exy_[pml_idx] 
          + com.get_e_y_coef2(jt) * grid.Cby_[mid] 
          * (grid.hz_[grid.pi(it, jt, kt)] 
             - grid.hz_[grid.pi(it, jt-1, kt)]);
        
        grid.ex_[grid_idx] = exz_[pml_idx] + exy_[pml_idx];
      }
}

void Pml::pml_update_ey(Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  region_t pml_r; // Region in the PML to update

  pml_r.xmin = pml_r.ymin = pml_r.zmin = 0;

  pml_r.xmax = grid_ey_r_.xmax - grid_ey_r_.xmin;
  pml_r.ymax = grid_ey_r_.ymax - grid_ey_r_.ymin;
  pml_r.zmax = grid_ey_r_.zmax - grid_ey_r_.zmin;

  PmlCommon &com = grid.get_pml_common();

  for(i = pml_r.xmin, it = grid_ey_r_.xmin; it < grid_ey_r_.xmax; i++, it++)
    for(j = pml_r.ymin, jt = grid_ey_r_.ymin; jt < grid_ey_r_.ymax; j++, jt++)
      for(k = pml_r.zmin, kt = grid_ey_r_.zmin; kt < grid_ey_r_.zmax; k++, kt++)
      {
        grid_idx = grid.pi(it, jt, kt);
        pml_idx = pi(i, j, k);

        mid = grid.material_[grid_idx];

        eyx_[pml_idx] = 
          com.get_e_x_coef1(it) * grid.Ca_[mid] * eyx_[pml_idx] 
          + com.get_e_x_coef2(it) * grid.Cbx_[mid] 
            * (grid.hz_[grid.pi(it-1, jt, kt)] 
               - grid.hz_[grid.pi(it, jt, kt)]);
        
        eyz_[pml_idx] = 
          com.get_e_z_coef1(kt) * grid.Ca_[mid] * eyz_[pml_idx] 
          + com.get_e_z_coef2(kt) * grid.Cbz_[mid] 
            * (grid.hx_[grid.pi(it, jt, kt)] 
               - grid.hx_[grid.pi(it, jt, kt-1)]);
        
        grid.ey_[grid_idx] = eyx_[pml_idx] + eyz_[pml_idx];
      }
}

void Pml::pml_update_ez(Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  region_t pml_r; // Region in the PML to update

  pml_r.xmin = pml_r.ymin = pml_r.zmin = 0;

  pml_r.xmax = grid_ez_r_.xmax - grid_ez_r_.xmin;
  pml_r.ymax = grid_ez_r_.ymax - grid_ez_r_.ymin;
  pml_r.zmax = grid_ez_r_.zmax - grid_ez_r_.zmin;

  PmlCommon &com = grid.get_pml_common();

  for(i = pml_r.xmin, it = grid_ez_r_.xmin; it < grid_ez_r_.xmax; i++, it++)
    for(j = pml_r.ymin, jt = grid_ez_r_.ymin; jt < grid_ez_r_.ymax; j++, jt++)
      for(k = pml_r.zmin, kt = grid_ez_r_.zmin; kt < grid_ez_r_.zmax; k++, kt++)
      {
        grid_idx = grid.pi(it, jt, kt);
        pml_idx = pi(i, j, k);

        mid = grid.material_[grid_idx];

        ezy_[pml_idx] = 
          com.get_e_y_coef1(jt) * grid.Ca_[mid] * ezy_[pml_idx] 
          + com.get_e_y_coef2(jt) * grid.Cby_[mid] 
          * (grid.hx_[grid.pi(it, jt-1, kt)] 
             - grid.hx_[grid.pi(it, jt, kt)]);
        
        ezx_[pml_idx] = 
          com.get_e_x_coef1(it) * grid.Ca_[mid] * ezx_[pml_idx] 
          + com.get_e_x_coef2(it) * grid.Cbx_[mid] 
          * (grid.hy_[grid.pi(it, jt, kt)] 
             - grid.hy_[grid.pi(it-1, jt, kt)]);
        
        grid.ez_[grid_idx] = ezx_[pml_idx] + ezy_[pml_idx];
      }
}

void Pml::pml_update_hx(Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  region_t pml_r; // Region in the PML to update

  pml_r.xmin = pml_r.ymin = pml_r.zmin = 0;

  pml_r.xmax = grid_hx_r_.xmax - grid_hx_r_.xmin;
  pml_r.ymax = grid_hx_r_.ymax - grid_hx_r_.ymin;
  pml_r.zmax = grid_hx_r_.zmax - grid_hx_r_.zmin;

  PmlCommon &com = grid.get_pml_common();

  for(i = pml_r.xmin, it = grid_hx_r_.xmin; it < grid_hx_r_.xmax; i++, it++)
    for(j = pml_r.ymin, jt = grid_hx_r_.ymin; jt < grid_hx_r_.ymax; j++, jt++)
      for(k = pml_r.zmin, kt = grid_hx_r_.zmin; kt < grid_hx_r_.zmax; k++, kt++)
      {
        grid_idx = grid.pi(it, jt, kt);
        pml_idx = pi(i, j, k);

        mid = grid.material_[grid_idx];

        hxz_[pml_idx] = 
          com.get_h_z_coef1(kt) * grid.Da_[mid] * hxz_[pml_idx] 
          + com.get_h_z_coef2(kt) * grid.Dbz_[mid] 
          * (grid.ey_[grid.pi(it, jt, kt+1)] 
             - grid.ey_[grid.pi(it, jt, kt)]);
        
        hxy_[pml_idx] = 
          com.get_h_y_coef1(jt) * grid.Da_[mid] * hxy_[pml_idx] 
          + com.get_h_y_coef2(jt) * grid.Dby_[mid] 
          * (grid.ez_[grid.pi(it, jt, kt)] 
             - grid.ez_[grid.pi(it, jt+1, kt)]);
        
        grid.hx_[grid_idx] = hxz_[pml_idx] + hxy_[pml_idx];
      }
}

void Pml::pml_update_hy(Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  region_t pml_r; // Region in the PML to update

  pml_r.xmin = pml_r.ymin = pml_r.zmin = 0;

  pml_r.xmax = grid_hy_r_.xmax - grid_hy_r_.xmin;
  pml_r.ymax = grid_hy_r_.ymax - grid_hy_r_.ymin;
  pml_r.zmax = grid_hy_r_.zmax - grid_hy_r_.zmin;

  PmlCommon &com = grid.get_pml_common();

  for(i = pml_r.xmin, it = grid_hy_r_.xmin; it < grid_hy_r_.xmax; i++, it++)
    for(j = pml_r.ymin, jt = grid_hy_r_.ymin; jt < grid_hy_r_.ymax; j++, jt++)
      for(k = pml_r.zmin, kt = grid_hy_r_.zmin; kt < grid_hy_r_.zmax; k++, kt++)
      {
        grid_idx = grid.pi(it, jt, kt);
        pml_idx = pi(i, j, k);

        mid = grid.material_[grid_idx];

        hyx_[pml_idx] = 
          com.get_h_x_coef1(it) * grid.Da_[mid] * hyx_[pml_idx] 
          + com.get_h_x_coef2(it) * grid.Dbx_[mid] 
          * (grid.ez_[grid.pi(it+1, jt, kt)] 
             - grid.ez_[grid.pi(it, jt, kt)]);
        
        hyz_[pml_idx] = 
          com.get_h_z_coef1(kt) * grid.Da_[mid] * hyz_[pml_idx] 
          + com.get_h_z_coef2(kt) * grid.Dbz_[mid] 
          * (grid.ex_[grid.pi(it, jt, kt)] 
             - grid.ex_[grid.pi(it, jt, kt+1)]);
        
        grid.hy_[grid_idx] = hyx_[pml_idx] + hyz_[pml_idx];
      }
}

void Pml::pml_update_hz(Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  region_t pml_r; // Region in the PML to update

  pml_r.xmin = pml_r.ymin = pml_r.zmin = 0;

  pml_r.xmax = grid_hz_r_.xmax - grid_hz_r_.xmin;
  pml_r.ymax = grid_hz_r_.ymax - grid_hz_r_.ymin;
  pml_r.zmax = grid_hz_r_.zmax - grid_hz_r_.zmin;

  PmlCommon &com = grid.get_pml_common();

  for(i = pml_r.xmin, it = grid_hz_r_.xmin; it < grid_hz_r_.xmax; i++, it++)
    for(j = pml_r.ymin, jt = grid_hz_r_.ymin; jt < grid_hz_r_.ymax; j++, jt++)
      for(k = pml_r.zmin, kt = grid_hz_r_.zmin; kt < grid_hz_r_.zmax; k++, kt++)
      {
        grid_idx = grid.pi(it, jt, kt);
        pml_idx = pi(i, j, k);

        mid = grid.material_[grid_idx];

        hzy_[pml_idx] = 
          com.get_h_y_coef1(jt) * grid.Da_[mid] * hzy_[pml_idx] 
          + com.get_h_y_coef2(jt) * grid.Dby_[mid] 
          * (grid.ex_[grid.pi(it, jt+1, kt)] 
             - grid.ex_[grid.pi(it, jt, kt)]);
        
        hzx_[pml_idx] = 
          com.get_h_x_coef1(it) * grid.Da_[mid] * hzx_[pml_idx] 
          + com.get_h_x_coef2(it) * grid.Dbx_[mid] 
          * (grid.ey_[grid.pi(it, jt, kt)] 
             - grid.ey_[grid.pi(it+1, jt, kt)]);
        
        grid.hz_[grid_idx] = hzx_[pml_idx] + hzy_[pml_idx];
      }
}

float Pml::sigma_over_eps_int(float x)
{
  if (!geometric_profile_)
  {
    if (x <= 0.0)
      return 0.0;

    else if (x <= delta_bndy_)
      return ratio_m_ * delta_bndy_ / (exponent_n_ + 1.0) 
        * (1.0 - pow(static_cast<float>((delta_bndy_ - x) / delta_bndy_), 
                     static_cast<float>(exponent_n_ + 1.0)));
    
    else
      return ratio_m_ * delta_bndy_ / (exponent_n_ + 1.0);
  }
  else
  {
    if (x <= 0.0)
      return 0.0;

    else if (x <= delta_bndy_)
      return ratio_m_ * delta_bndy_ * pow(g_, delta_bndy_ / geometric_delta_)
        / log(g_) * (1.0 - pow(g_, -x / geometric_delta_));

    else 
      return ratio_m_ * delta_bndy_ 
        * (pow(g_, delta_bndy_ / geometric_delta_) - 1.0)
        / log(g_);
  }
}

void Pml::add_sd_bcs(SubdomainBc *sd, Face pmlface, Face sdface)
{
  RxTxData ret;
  unsigned int s_idx = 0, r_idx = 0;

  switch (sdface) {
  case BACK:
    s_idx = pi(1, 0, 0);
    r_idx = pi(0, 0, 0);
    ret.set_datatype(xz_plane_);
    break;
  case FRONT:
    s_idx = pi(pml_r_.xmax - 2, 0, 0);
    r_idx = pi(pml_r_.xmax - 1, 0, 0);
    ret.set_datatype(xz_plane_);
    break;
  case LEFT:
    s_idx = pi(0, 1, 0);
    r_idx = pi(0, 0, 0);
    ret.set_datatype(yz_plane_);
    break;
  case RIGHT:
    s_idx = pi(0, pml_r_.ymax - 2, 0);
    r_idx = pi(0, pml_r_.ymax - 1, 0);
    ret.set_datatype(yz_plane_);
    break;
  case BOTTOM:
    s_idx = pi(0, 0, 0);
    r_idx = pi(0, 0, 1);
    ret.set_datatype(xy_plane_);
    break;
  case TOP:
    s_idx = pi(0, 0, pml_r_.zmax - 2);
    r_idx = pi(0, 0, pml_r_.zmax - 1);
    ret.set_datatype(xy_plane_);
    break;
  }

  ret.set_field_type(E);

  ret.set_rx_ptr(&(exy_[r_idx]));
  ret.set_tx_ptr(&(exy_[s_idx]));
  sd->add_tx_rx_data(ret);

  ret.set_rx_ptr(&(exz_[r_idx]));
  ret.set_tx_ptr(&(exz_[s_idx]));
  sd->add_tx_rx_data(ret);

  ret.set_rx_ptr(&(eyx_[r_idx]));
  ret.set_tx_ptr(&(eyx_[s_idx]));
  sd->add_tx_rx_data(ret);

  ret.set_rx_ptr(&(eyz_[r_idx]));
  ret.set_tx_ptr(&(eyz_[s_idx]));
  sd->add_tx_rx_data(ret);

  ret.set_rx_ptr(&(ezx_[r_idx]));
  ret.set_tx_ptr(&(ezx_[s_idx]));
  sd->add_tx_rx_data(ret);

  ret.set_rx_ptr(&(ezy_[r_idx]));
  ret.set_tx_ptr(&(ezy_[s_idx]));
  sd->add_tx_rx_data(ret);

  ret.set_field_type(H);

  ret.set_rx_ptr(&(hxy_[r_idx]));
  ret.set_tx_ptr(&(hxy_[s_idx]));
  sd->add_tx_rx_data(ret);

  ret.set_rx_ptr(&(hxz_[r_idx]));
  ret.set_tx_ptr(&(hxz_[s_idx]));
  sd->add_tx_rx_data(ret);

  ret.set_rx_ptr(&(hyx_[r_idx]));
  ret.set_tx_ptr(&(hyx_[s_idx]));
  sd->add_tx_rx_data(ret);

  ret.set_rx_ptr(&(hyz_[r_idx]));
  ret.set_tx_ptr(&(hyz_[s_idx]));
  sd->add_tx_rx_data(ret);

  ret.set_rx_ptr(&(hzx_[r_idx]));
  ret.set_tx_ptr(&(hzx_[s_idx]));
  sd->add_tx_rx_data(ret);

  ret.set_rx_ptr(&(hzy_[r_idx]));
  ret.set_tx_ptr(&(hzy_[s_idx]));
  sd->add_tx_rx_data(ret);

}

BoundaryCondition Pml::get_type()
{
  return PML;
}
