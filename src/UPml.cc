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
  : d_(0), h_(0), sigmas_(0), poly_order_(4), C1_(0), C2_(0)
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

    switch (face)
    {
    case FRONT:
    case BACK:
      if (grid.get_boundary(LEFT).get_type() == UPML)
      {
        grid_ey_r_.ymin += grid.get_boundary(LEFT).get_thickness();
        grid_hy_r_.ymin += grid.get_boundary(LEFT).get_thickness() + 1;
      }

      if (grid.get_boundary(RIGHT).get_type() == UPML)
      {
        grid_ey_r_.ymax -= grid.get_boundary(LEFT).get_thickness();
        grid_hy_r_.ymax -= grid.get_boundary(LEFT).get_thickness() + 1;
      }

      if (grid.get_boundary(TOP).get_type() == UPML)
      {
        grid_ez_r_.zmax -= grid.get_boundary(TOP).get_thickness();
        grid_hz_r_.zmax -= grid.get_boundary(TOP).get_thickness() + 1;
      }

      if (grid.get_boundary(BOTTOM).get_type() == UPML)
      {
        grid_ez_r_.zmin += grid.get_boundary(BOTTOM).get_thickness();
        grid_hz_r_.zmin += grid.get_boundary(BOTTOM).get_thickness() + 1;
      }
      break;

    case LEFT:
    case RIGHT:
      if (grid.get_boundary(FRONT).get_type() == UPML)
      {
        unsigned int thickness = grid.get_boundary(FRONT).get_thickness();
        grid_ex_r_.xmax -= thickness;
        grid_hx_r_.xmax -= thickness + 1;
        grid_ez_r_.xmax -= thickness;
        grid_hz_r_.xmax -= thickness + 1;
      }

      if (grid.get_boundary(BACK).get_type() == UPML)
      {
        unsigned int thickness = grid.get_boundary(BACK).get_thickness();
        grid_ex_r_.xmin += thickness;
        grid_hx_r_.xmin += thickness + 1;
        grid_ez_r_.xmin += thickness;
        grid_hz_r_.xmin += thickness + 1;
      }

      if (grid.get_boundary(TOP).get_type() == UPML)
      {
        unsigned int thickness = grid.get_boundary(TOP).get_thickness();
        grid_ez_r_.zmax -= thickness;
        grid_hz_r_.zmax -= thickness + 1;
        grid_ex_r_.zmax -= thickness;
        grid_hx_r_.zmax -= thickness + 1;
      }

      if (grid.get_boundary(BOTTOM).get_type() == UPML)
      {
        unsigned int thickness = grid.get_boundary(BOTTOM).get_thickness();
        grid_ez_r_.zmin += thickness;
        grid_hz_r_.zmin += thickness + 1;
        grid_ex_r_.zmin += thickness;
        grid_hx_r_.zmin += thickness + 1;
      }
      break;

    case TOP:
    case BOTTOM:
      if (grid.get_boundary(LEFT).get_type() == UPML)
      {
        unsigned int thickness = grid.get_boundary(LEFT).get_thickness();
        grid_ey_r_.ymin += thickness;
        grid_hy_r_.ymin += thickness + 1;
        grid_ex_r_.ymin += thickness;
        grid_hx_r_.ymin += thickness + 1;
      }

      if (grid.get_boundary(RIGHT).get_type() == UPML)
      {
        unsigned int thickness = grid.get_boundary(LEFT).get_thickness();
        grid_ey_r_.ymax -= thickness;
        grid_hy_r_.ymax -= thickness + 1;
        grid_ex_r_.ymax -= thickness;
        grid_hx_r_.ymax -= thickness + 1;
      }

      if (grid.get_boundary(FRONT).get_type() == UPML)
      {
        unsigned int thickness = grid.get_boundary(FRONT).get_thickness();
        grid_ex_r_.xmax -= thickness;
        grid_hx_r_.xmax -= thickness + 1;
        grid_ey_r_.xmax -= thickness;
        grid_hy_r_.xmax -= thickness;
      }

      if (grid.get_boundary(BACK).get_type() == UPML)
      {
        unsigned int thickness = grid.get_boundary(BACK).get_thickness();
        grid_ex_r_.xmin += thickness;
        grid_hx_r_.xmin += thickness + 1;
        grid_ey_r_.xmin += thickness;
        grid_hy_r_.xmin += thickness;
      }
      break;
    }
  }
}

void UPml::init(const Grid &grid, Face face)
{
  d_file.open("upml_d.txt");
  h_file.open("upml_h.txt");

  compute_regions(face, grid);

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

  sigmas_ = new mat_coef_t[thickness_ + 1];

  const MaterialLib *mlib = &(grid.get_material_lib());
  
  vector<Material>::const_iterator iter = mlib->get_material_iter_begin();
  vector<Material>::const_iterator iter_e = mlib->get_material_iter_end();
  int index = 0;
  
  ++index;
  ++iter;

  while (iter != iter_e) 
  {
    mat_prop_t eps = (*iter).get_epsilon() * EPS_0;
    mat_prop_t sig = (*iter).get_sigma();
    mat_prop_t mu = (*iter).get_mu() * MU_0;
    //mat_prop_t sigs = (*iter).get_sigma_star();

    mat_coef_t sigma_max = (poly_order_ + 1) 
      / (150 * PI * delta * sqrt((*iter).get_epsilon()));

    for (int idx = 0; idx < thickness_ + 1; idx++)
    {
      sigmas_[idx] = sigma_max * pow(static_cast<float>(idx + 1), 
                                     static_cast<float>(poly_order_)) 
        / pow(static_cast<float>(thickness_),
              static_cast<float>(poly_order_));
    }

//     cerr << "Sigmas: ";
//     for (int idx = 0; idx < thickness_ + 1; idx++)
//     {
//       cerr << sigmas_[idx] << " ";
//     }
//     cerr << endl;

    if (BACK == face || BOTTOM == face || LEFT == face)
    {
      for (int idx = 0; idx < (thickness_ + 1)/2; idx++)
      {
        mat_coef_t temp = sigmas_[idx];
        sigmas_[idx] = sigmas_[thickness_ - idx];
        sigmas_[thickness_ - idx] = temp;
      }
    }

//     cerr << "Sigmas: ";
//     for (int idx = 0; idx < thickness_ + 1; idx++)
//     {
//       cerr << sigmas_[idx] << " ";
//     }
//     cerr << endl;


    for (int idx = 0; idx < thickness_ + 1; idx++)
    {
      sig = sigmas_[idx];

//       C1_[idx][index] = (1 + (sig * grid.get_deltat())/(2*EPS_0)) / eps;
//       C2_[idx][index] = (1 - (sig * grid.get_deltat())/(2*EPS_0)) / eps;

//       D1_[idx][index] = (1 + (sig * grid.get_deltat())/(2*EPS_0)) / mu;
//       D2_[idx][index] = (1 - (sig * grid.get_deltat())/(2*EPS_0)) / mu;

      C1_[idx][index] = (1 + (sig * grid.get_deltat() * 0.5)/eps);
      C2_[idx][index] = (1 - (sig * grid.get_deltat() * 0.5)/eps);

      D1_[idx][index] = (1 + (sig * grid.get_deltat() * 0.5)/mu);
      D2_[idx][index] = (1 - (sig * grid.get_deltat() * 0.5)/mu);

//       cerr << "idx: " << idx << ", index: " << index 
//            << ", C1 = " << C1_[idx][index]
//            << ", C2 = " << C2_[idx][index]
//            << ", D1 = " << D1_[idx][index]
//            << ", D2 = " << D2_[idx][index] << endl;
    }

    ++iter;
    ++index;      
  }  

//   cout << "UPML Update region for face " << face << ":"
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

}

void UPml::alloc_coefs(unsigned int num_materials)
{
  free_coefs();

  C1_ = new mat_coef_t*[thickness_ + 1];
  C2_ = new mat_coef_t*[thickness_ + 1];

  D1_ = new mat_coef_t*[thickness_ + 1];
  D2_ = new mat_coef_t*[thickness_ + 1];

  if (!C1_ || !C2_ || !D1_ || !D2_)
  {
    free_coefs();
    throw MemoryException();
  }

  for (int idx = 0; idx <= thickness_; idx++)
  {
    C1_[idx] = new mat_coef_t[num_materials];
    C2_[idx] = new mat_coef_t[num_materials];

    D1_[idx] = new mat_coef_t[num_materials];
    D2_[idx] = new mat_coef_t[num_materials];

    if (!C1_[idx] || !C2_[idx] || !D1_[idx] || !D2_[idx])
    {
      free_coefs();
      throw MemoryException();
    }
  }
}

void UPml::free_coefs()
{
  if (C1_)
  {
    for (int idx = 0; idx <= thickness_; idx++)
    {
      if (C1_[idx])
        delete[] C1_[idx];
      if (C2_[idx])
        delete[] C2_[idx];

      if (D1_[idx])
        delete[] D1_[idx];
      if (D2_[idx])
        delete[] D2_[idx];
    }
    delete[] C1_;
    C1_ = 0;
    delete[] C2_;

    delete[] D1_;
    delete[] D2_;
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
    for(sig_idx = 0, i = pml_r.xmin, it = grid_ex_r_.xmin; 
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
            d_temp = d_[pml_idx] + grid.get_deltat() 
              * ((grid.hz_[grid.pi(it, jt, kt)] 
                  - grid.hz_[grid.pi(it, jt-1, kt)]) / grid.get_deltay()
                 + (grid.hy_[grid_idx - 1] 
                    - grid.hy_[grid_idx]) / grid.get_deltaz());
            
            grid.ex_[grid_idx] = grid.ex_[grid_idx]
              + ( C1_[sig_idx][mid] * d_temp
                  - C2_[sig_idx][mid] * d_[pml_idx]);
           
            if (i == 2 && j ==9 && k == 9)
              d_file << "ex_ = " << grid.ex_[grid_idx] << ", d_ = " 
                     << d_temp << "\n\thz1_ = " << grid.hz_[grid_idx]
                     << " hz2_ = " << grid.hz_[grid.pi(it, jt-1, kt)] 
                     << "\n\thy1_ = " << grid.hy_[grid_idx - 1]
                     << " hy2_ = " << grid.hy_[grid_idx] << endl;

            d_[pml_idx] = d_temp;
            
          } else {
            grid.ex_[grid_idx] = grid.Ca_[mid] * grid.ex_[grid_idx]
              + grid.Cby_[mid] 
              * (grid.hz_[grid.pi(it, jt, kt)] 
                 - grid.hz_[grid.pi(it, jt-1, kt)])
              + grid.Cbz_[mid] 
              * (grid.hy_[grid_idx - 1] 
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
  unsigned int grid_idx, pml_idx, mid, sig_idx = 0; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_ey_r_); 
  
  field_t d_temp = 0;

#ifdef USE_OPENMP
#pragam omp parallel private(mid, grid_idx, pml_idx, sig_idx, i, j, k, it, jt, kt, d_temp)
  {
#pragam omp for
#endif
    for(i = pml_r.xmin, it = grid_ey_r_.xmin; 
        it < grid_ey_r_.xmax; i++, it++)
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
          
          if (pml) 
          {
            d_temp = d_[pml_idx] + grid.get_deltat() 
              * ((grid.hx_[grid_idx] 
                  - grid.hx_[grid.pi(it, jt, kt-1)]) / grid.get_deltaz()
                 + (grid.hz_[grid.pi(it-1, jt, kt)] 
                    - grid.hz_[grid_idx]) / grid.get_deltax());
            
            grid.ey_[grid_idx] = grid.ey_[grid_idx]
              + ( C1_[sig_idx][mid] * d_temp
                  - C2_[sig_idx][mid] * d_[pml_idx]);
                                     
            d_[pml_idx] = d_temp;
                                     
          } else {
            grid.ey_[grid_idx] = grid.Ca_[mid] * grid.ey_[grid_idx]
              + grid.Cbz_[mid] * (grid.hx_[grid_idx] 
                                  - grid.hx_[grid.pi(it, jt, kt-1)])
              + grid.Cbx_[mid] * (grid.hz_[grid.pi(it-1, jt, kt)] 
                                  - grid.hz_[grid_idx]);
          }
        }
      }
    }
#ifdef USE_OPENMP
  }
#endif
}

void UPml::update_ez(Grid &grid, bool pml) 
{
  unsigned int grid_idx, pml_idx, mid, sig_idx = 0; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_ez_r_); 
  
  field_t d_temp = 0;

#ifdef USE_OPENMP
#pragam omp parallel private(mid, grid_idx, pml_idx, sig_idx, i, j, k, it, jt, kt, d_temp)
  {
#pragam omp for
#endif
    for(i = pml_r.xmin, it = grid_ez_r_.xmin; 
        it < grid_ez_r_.xmax; i++, it++)
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
          
          if (pml) 
          {
            d_temp = d_[pml_idx] + grid.get_deltat() 
              * ((grid.hy_[grid.pi(it, jt, kt)] 
                  - grid.hy_[grid.pi(it-1, jt, kt)]) / grid.get_deltax()
                 + (grid.hx_[grid.pi(it, jt - 1, kt)] 
                    - grid.hx_[grid_idx]) / grid.get_deltay());
            
            grid.ez_[grid_idx] = grid.ez_[grid_idx]
              + ( C1_[sig_idx][mid] * d_temp
                  - C2_[sig_idx][mid] * d_[pml_idx]);
                                     
            d_[pml_idx] = d_temp;
                                     
          } else {
            grid.ez_[grid_idx] = grid.Ca_[mid] * grid.ez_[grid_idx]
              + grid.Cbx_[mid] * (grid.hy_[grid.pi(it, jt, kt)] 
                                  - grid.hy_[grid.pi(it-1, jt, kt)])
              + grid.Cby_[mid] * (grid.hx_[grid.pi(it, jt - 1, kt)] 
                                  - grid.hx_[grid_idx]);
          }
        }
      }
    }
#ifdef USE_OPENMP
  }
#endif
}

void UPml::update_hx(Grid &grid, bool pml)
{
  unsigned int grid_idx, pml_idx, mid, sig_idx = 0; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_hx_r_); 
  
  field_t h_temp = 0;

#ifdef USE_OPENMP
#pragam omp parallel private(mid, grid_idx, pml_idx, sig_idx, i, j, k, it, jt, kt, h_temp)
  {
#pragam omp for
#endif
    for(sig_idx = 0, i = pml_r.xmin, it = grid_hx_r_.xmin; 
        it < grid_hx_r_.xmax; i++, it++, sig_idx++)
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
          
          if (pml) 
          {
            h_temp = h_[pml_idx] - grid.get_deltat() 
              * ((grid.ez_[grid.pi(it, jt, kt)] 
                  - grid.ez_[grid.pi(it, jt+1, kt)]) / grid.get_deltay()
                 + (grid.ey_[grid_idx + 1] 
                    - grid.ey_[grid_idx]) / grid.get_deltaz());
            
            grid.hx_[grid_idx] = grid.hx_[grid_idx]
              + ( D1_[sig_idx][mid] * h_temp
                  - D2_[sig_idx][mid] * h_[pml_idx]);

            if (i == 2 && j ==9 && k == 9)
              h_file << "hx_ = " << grid.hx_[grid_idx] << ", h_ = " 
                     << h_temp << "\n\tez1_ = " << grid.ez_[grid_idx]
                     << " ez2_ = " << grid.ez_[grid.pi(it, jt+1, kt)] 
                     << "\n\tey1_ = " << grid.ey_[grid_idx + 1]
                     << " ey2_ = " << grid.ey_[grid_idx] << endl;

            h_[pml_idx] = h_temp;
                                     
          } else {
            grid.hx_[grid_idx] = grid.Da_[mid] * grid.hx_[grid_idx]
              + grid.Dby_[mid] * (grid.ez_[grid.pi(it, jt, kt)] 
                             - grid.ez_[grid.pi(it, jt+1, kt)])
              + grid.Dbz_[mid] * (grid.ey_[grid_idx + 1] 
                             - grid.ey_[grid_idx]);
          }
        }
      }
    }
#ifdef USE_OPENMP
  }
#endif
}

void UPml::update_hy(Grid &grid, bool pml)
{
  unsigned int grid_idx, pml_idx, mid, sig_idx = 0; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_hy_r_); 
  
  field_t h_temp = 0;

#ifdef USE_OPENMP
#pragam omp parallel private(mid, grid_idx, pml_idx, sig_idx, i, j, k, it, jt, kt, h_temp)
  {
#pragam omp for
#endif
    for(i = pml_r.xmin, it = grid_hy_r_.xmin; 
        it < grid_hy_r_.xmax; i++, it++)
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
          
          if (pml) 
          {
            h_temp = h_[pml_idx] - grid.get_deltat() 
              * ((grid.ex_[grid_idx] 
                  - grid.ex_[grid.pi(it, jt, kt+1)]) / grid.get_deltaz()
                 + (grid.ez_[grid.pi(it+1, jt, kt)] 
                    - grid.ez_[grid_idx]) / grid.get_deltax());
            
            grid.hy_[grid_idx] = grid.hy_[grid_idx]
              + ( D1_[sig_idx][mid] * h_temp
                  - D2_[sig_idx][mid] * h_[pml_idx]);
                                     
            h_[pml_idx] = h_temp;
                                     
          } else {
            grid.hy_[grid_idx] = grid.Da_[mid] * grid.hy_[grid_idx]
              + grid.Dbz_[mid] * (grid.ex_[grid_idx] 
                                  - grid.ex_[grid.pi(it, jt, kt+1)])
              + grid.Dbx_[mid] * (grid.ez_[grid.pi(it+1, jt, kt)] 
                                  - grid.ez_[grid_idx]);
          }
        }
      }
    }
#ifdef USE_OPENMP
  }
#endif
}

void UPml::update_hz(Grid &grid, bool pml)
{
  unsigned int grid_idx, pml_idx, mid, sig_idx = 0; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  // Region in the PML to update
  region_t pml_r = find_local_region(grid_hz_r_); 
  
  field_t h_temp = 0;

#ifdef USE_OPENMP
#pragam omp parallel private(mid, grid_idx, pml_idx, sig_idx, i, j, k, it, jt, kt, h_temp)
  {
#pragam omp for
#endif
    for(i = pml_r.xmin, it = grid_hz_r_.xmin; 
        it < grid_hz_r_.xmax; i++, it++)
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
          
          if (pml) 
          {
            h_temp = h_[pml_idx] - grid.get_deltat() 
              * ((grid.ey_[grid.pi(it, jt, kt)] 
                  - grid.ey_[grid.pi(it+1, jt, kt)]) / grid.get_deltax()
                 + (grid.ex_[grid.pi(it, jt + 1, kt)] 
                    - grid.ex_[grid_idx]) / grid.get_deltay());
            
            grid.hz_[grid_idx] = grid.hz_[grid_idx]
              + ( D1_[sig_idx][mid] * h_temp
                  - D2_[sig_idx][mid] * h_[pml_idx]);
                                     
            h_[pml_idx] = h_temp;
                                     
          } else {
            grid.hz_[grid_idx] = grid.Da_[mid] * grid.hz_[grid_idx]
              + grid.Dbx_[mid] * (grid.ey_[grid.pi(it, jt, kt)] 
                                  - grid.ey_[grid.pi(it+1, jt, kt)])
              + grid.Dby_[mid] * (grid.ex_[grid.pi(it, jt + 1, kt)] 
                                  - grid.ex_[grid_idx]);
          }
        }
      }
    }
#ifdef USE_OPENMP
  }
#endif
}
