#include "Pml.hh"
#include "Grid.hh"
#include <exception>

using namespace std;

Pml::Pml()
  : exy_(0), exz_(0), eyx_(0), eyz_(0), ezx_(0), ezy_(0),
  hxy_(0), hxz_(0), hyx_(0), hyz_(0), hzx_(0), hzy_(0), alloced_(false)
{}

Pml::~Pml()
{
  free_pml_fields();
}

void Pml::alloc_pml_fields(Face face, Grid &grid)
{
  if (alloced_)
    return; 

  if (thickness_ == 0)
    throw exception();
  
  region_t r = find_face(face, grid);
  pml_r_ = r;

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
    throw exception(); // Insufficent memory
  }

}

void Pml::set_thickness(unsigned int thickness)
{
  thickness_ = thickness;
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
}

void Pml::apply(Face face, Grid &grid)
{
  region_t grid_r = find_face(face, grid);
  region_t e_grid_r = grid_r; 

  // Modify the grid region so that the outer walls are not computed;
  // that they be electric walls. 
  switch (face)
  {
  case FRONT:
    e_grid_r.xmax--;
    break;
  case BACK:
    e_grid_r.xmin++;
    break;

  case LEFT:
    e_grid_r.ymin++;
    break;

  case RIGHT:
    e_grid_r.ymax--;
    break;

  case TOP:
    e_grid_r.zmax--;
    break;

  case BOTTOM:
    e_grid_r.zmin++;
    break;
  }

  pml_update_hx(grid_r, grid);
  pml_update_hy(grid_r, grid);
  pml_update_hz(grid_r, grid);

  pml_update_ex(e_grid_r, grid);
  pml_update_ey(e_grid_r, grid);
  pml_update_ez(e_grid_r, grid);
}

void Pml::pml_update_ex(const region_t &grid_r, Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */
  region_t tmp_r;

  PmlCommon &com = grid.get_pml_common();

  for(i = pml_r_.xmin, it = grid_r.xmin; it <= grid_r.xmax; i++, it++)
    for(j = pml_r_.ymin, jt = grid_r.ymin; jt <= grid_r.ymax-1; j++, jt++)
      for(k = pml_r_.zmin, kt = grid_r.zmin; kt <= grid_r.zmax-1; k++, kt++)
      {
        grid_idx = grid.pi(it, jt, kt);
        pml_idx = pi(i, j, k);

        mid = grid.material_[grid_idx];

        hxz_[pml_idx] = com.h_z_coef1_[kt] * grid.Da_[mid] *
          hxz_[pml_idx] + com.h_z_coef2_[kt] *
          grid.Dbz_[mid] * (grid.ey_[grid.pi(it, jt, kt+1)] 
                           - grid.ey_[grid.pi(it, jt, kt)]);
        
        hxy_[pml_idx] = com.h_y_coef1_[jt] * grid.Da_[mid] *
          hxy_[pml_idx] - com.h_y_coef2_[jt] *
          grid.Dby_[mid] * (grid.ez_[pi(it, jt+1, kt)] 
                           - grid.ez_[pi(it, jt, kt)]);
        
        grid.hx_[grid_idx] = hxz_[pml_idx] + hxy_[pml_idx];
      }
}

void Pml::pml_update_ey(const region_t &grid_r, Grid &grid)
{

}

void Pml::pml_update_ez(const region_t &grid_r, Grid &grid)
{

}

void Pml::pml_update_hx(const region_t &grid_r, Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  PmlCommon &com = grid.get_pml_common();

  for(i = pml_r_.xmin, it = grid_r.xmin; it <= grid_r.xmax; i++, it++)
    for(j = pml_r_.ymin, jt = grid_r.ymin; jt <= grid_r.ymax-1; j++, jt++)
      for(k = pml_r_.zmin, kt = grid_r.zmin; kt <= grid_r.zmax-1; k++, kt++)
      {
        grid_idx = grid.pi(it, jt, kt);
        pml_idx = pi(i, j, k);

        mid = grid.material_[grid_idx];

        hxz_[pml_idx] = com.h_z_coef1_[kt] * grid.Da_[mid] *
          hxz_[pml_idx] + com.h_z_coef2_[kt] *
          grid.Dbz_[mid] * (grid.ey_[grid.pi(it, jt, kt+1)] 
                           - grid.ey_[grid.pi(it, jt, kt)]);
        
        hxy_[pml_idx] = com.h_y_coef1_[jt] * grid.Da_[mid] *
          hxy_[pml_idx] - com.h_y_coef2_[jt] *
          grid.Dby_[mid] * (grid.ez_[pi(it, jt+1, kt)] 
                           - grid.ez_[pi(it, jt, kt)]);
        
        grid.hx_[grid_idx] = hxz_[pml_idx] + hxy_[pml_idx];
      }
}

void Pml::pml_update_hy(const region_t &grid_r, Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  PmlCommon &com = grid.get_pml_common();

  for(i = pml_r_.xmin, it = grid_r.xmin; it <= grid_r.xmax - 1; i++, it++)
    for(j = pml_r_.ymin, jt = grid_r.ymin; jt <= grid_r.ymax; j++, jt++)
      for(k = pml_r_.zmin, kt = grid_r.zmin; kt <= grid_r.zmax-1; k++, kt++)
      {
        grid_idx = grid.pi(it, jt, kt);
        pml_idx = pi(i, j, k);

        mid = grid.material_[grid_idx];

        hyx_[pml_idx] = com.h_x_coef1_[it] * grid.Da_[mid] *
          hyx_[pml_idx] + com.h_x_coef2_[it] *
          grid.Dbx_[mid] * (grid.ez_[grid.pi(it+1, jt, kt)] 
                           - grid.ez_[grid.pi(it, jt, kt)]);
        
        hyz_[pml_idx] = com.h_z_coef1_[kt] * grid.Da_[mid] *
          hyz_[pml_idx] - com.h_z_coef2_[kt] *
          grid.Dbz_[mid] * (grid.ex_[pi(it, jt, kt+1)] 
                           - grid.ex_[pi(it, jt, kt)]);
        
        grid.hy_[grid_idx] = hyx_[pml_idx] + hyz_[pml_idx];
      }
}

void Pml::pml_update_hz(const region_t &grid_r, Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  PmlCommon &com = grid.get_pml_common();

  for(i = pml_r_.xmin, it = grid_r.xmin; it <= grid_r.xmax-1; i++, it++)
    for(j = pml_r_.ymin, jt = grid_r.ymin; jt <= grid_r.ymax-1; j++, jt++)
      for(k = pml_r_.zmin, kt = grid_r.zmin; kt <= grid_r.zmax; k++, kt++)
      {
        grid_idx = grid.pi(it, jt, kt);
        pml_idx = pi(i, j, k);

        mid = grid.material_[grid_idx];

        hzy_[pml_idx] = com.h_y_coef1_[jt] * grid.Da_[mid] *
          hzy_[pml_idx] + com.h_y_coef2_[jt] *
          grid.Dby_[mid] * (grid.ex_[grid.pi(it, jt+1, kt)] 
                           - grid.ex_[grid.pi(it, jt, kt)]);
        
        hzx_[pml_idx] = com.h_x_coef1_[it] * grid.Da_[mid] *
          hzx_[pml_idx] - com.h_x_coef2_[it] *
          grid.Dbx_[mid] * (grid.ey_[pi(it+1, jt, kt)] 
                           - grid.ey_[pi(it, jt, kt)]);
        
        grid.hz_[grid_idx] = hzx_[pml_idx] + hzy_[pml_idx];
      }
}
