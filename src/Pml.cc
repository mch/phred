#include "Pml.hh"
#include "Grid.hh"
#include "Exceptions.hh"

using namespace std;

Pml::Pml()
  : variation_(P), g_(0.0), nrml_refl_(1.0),
    ratio_m_(0.0), exponent_n_(0.0), delta_bndy_(0.0), 
    geometric_delta_(0.0), geometric_profile_(0),
    exy_(0), exz_(0), eyx_(0), eyz_(0), ezx_(0), ezy_(0),
    hxy_(0), hxz_(0), hyx_(0), hyz_(0), hzx_(0), hzy_(0), alloced_(false)
{}

Pml::Pml(PmlVariation variation, float g, float nrml_refl)
  : variation_(variation), g_(g), nrml_refl_(nrml_refl),
    ratio_m_(0.0), exponent_n_(0.0), delta_bndy_(0.0), 
    geometric_delta_(0.0), geometric_profile_(0),
    exy_(0), exz_(0), eyx_(0), eyz_(0), ezx_(0), ezy_(0),
    hxy_(0), hxz_(0), hyx_(0), hyz_(0), hzx_(0), hzy_(0), alloced_(false)
{}

Pml::Pml(char variation, float nrml_refl)
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
    throw MemoryException(); // Insufficent memory
  }

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
    ratio_m_ *= log(g_) / (pow(g_, static_cast<double>(thickness)) - 1.0);
    geometric_profile_ = 1;
    break;
  }

  cout << "Pml setup TEST results for face " << face << "\n" 
       << "\tratio_m_ = " << ratio_m_
       << "\n\texponent_n_ = " << exponent_n_
       << "\n\tdelta_bndy_ = " << delta_bndy_
       << "\n\tgeometric_delta_ = " << geometric_delta_
       << "\n\tgeomtric_profile_ = " << geometric_profile_ << endl;

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
  if (!alloced_)
    throw exception(); // PML must be set up before applying it. 

  region_t grid_r = find_face(face, grid);
  region_t e_grid_r = grid_r; 
  region_t e_pml_r = pml_r_;

  // Modify the grid region so that the outer walls are not computed;
  // that they be electric walls. WARNING! Have to change the PML region too!!
  switch (face)
  {
  case FRONT:
    e_grid_r.xmax--;
    e_pml_r.xmax--;
    break;
  case BACK:
    e_grid_r.xmin++;
    e_pml_r.xmin++;
    break;

  case LEFT:
    e_grid_r.ymin++;
    e_pml_r.ymin++;
    break;

  case RIGHT:
    e_grid_r.ymax--;
    e_pml_r.ymax--;
    break;

  case TOP:
    e_grid_r.zmax--;
    e_pml_r.zmax--;
    break;

  case BOTTOM:
    e_grid_r.zmin++;
    e_pml_r.zmin++;
    break;
  }

  pml_update_hx(grid_r, grid);
  pml_update_hy(grid_r, grid);
  pml_update_hz(grid_r, grid);

  pml_update_ex(e_pml_r, e_grid_r, grid);
  pml_update_ey(e_pml_r, e_grid_r, grid);
  pml_update_ez(e_pml_r, e_grid_r, grid);
}

void Pml::pml_update_ex(const region_t &pml_r, 
                        const region_t &grid_r, Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */
  region_t tmp_r;

  PmlCommon &com = grid.get_pml_common();

  for(i = pml_r.xmin, it = grid_r.xmin; it <= grid_r.xmax; i++, it++)
    for(j = pml_r.ymin+1, jt = grid_r.ymin+1; jt <= grid_r.ymax; j++, jt++)
      for(k = pml_r.zmin+1, kt = grid_r.zmin+1; kt <= grid_r.zmax; k++, kt++)
      {
        grid_idx = grid.pi(it, jt, kt);
        pml_idx = pi(i, j, k);

        mid = grid.material_[grid_idx];

        exz_[pml_idx] = com.e_z_coef1_[kt] * grid.Ca_[mid] *
          exz_[pml_idx] + com.e_z_coef2_[kt] *
          grid.Cbz_[mid] * (grid.hy_[grid.pi(it, jt, kt)] 
                           - grid.hy_[grid.pi(it, jt, kt-1)]);
        
        exy_[pml_idx] = com.e_y_coef1_[jt] * grid.Ca_[mid] *
          exy_[pml_idx] - com.e_y_coef2_[jt] *
          grid.Cby_[mid] * (grid.hz_[pi(it, jt, kt)] 
                           - grid.hz_[pi(it, jt-1, kt)]);
        
        grid.ex_[grid_idx] = exz_[pml_idx] + exy_[pml_idx];
      }
}

void Pml::pml_update_ey(const region_t &pml_r, 
                        const region_t &grid_r, Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  PmlCommon &com = grid.get_pml_common();

  for(i = pml_r.xmin+1, it = grid_r.xmin+1; it <= grid_r.xmax; i++, it++)
    for(j = pml_r.ymin, jt = grid_r.ymin; jt <= grid_r.ymax; j++, jt++)
      for(k = pml_r.zmin+1, kt = grid_r.zmin+1; kt <= grid_r.zmax; k++, kt++)
      {
        grid_idx = grid.pi(it, jt, kt);
        pml_idx = pi(i, j, k);

        mid = grid.material_[grid_idx];

        eyx_[pml_idx] = com.e_x_coef1_[it] * grid.Ca_[mid] *
          eyx_[pml_idx] + com.e_x_coef2_[it] *
          grid.Cbx_[mid] * (grid.hz_[grid.pi(it, jt, kt)] 
                           - grid.hz_[grid.pi(it-1, jt, kt)]);
        
        eyz_[pml_idx] = com.e_z_coef1_[kt] * grid.Ca_[mid] *
          eyz_[pml_idx] - com.e_z_coef2_[kt] *
          grid.Cbz_[mid] * (grid.hx_[pi(it, jt, kt)] 
                           - grid.hx_[pi(it, jt, kt-1)]);
        
        grid.ey_[grid_idx] = eyx_[pml_idx] + eyz_[pml_idx];
      }
}

void Pml::pml_update_ez(const region_t &pml_r, 
                        const region_t &grid_r, Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  PmlCommon &com = grid.get_pml_common();

  for(i = pml_r.xmin+1, it = grid_r.xmin+1; it <= grid_r.xmax; i++, it++)
    for(j = pml_r.ymin+1, jt = grid_r.ymin+1; jt <= grid_r.ymax; j++, jt++)
      for(k = pml_r.zmin, kt = grid_r.zmin; kt <= grid_r.zmax; k++, kt++)
      {
        grid_idx = grid.pi(it, jt, kt);
        pml_idx = pi(i, j, k);

        mid = grid.material_[grid_idx];

        ezy_[pml_idx] = com.e_y_coef1_[jt] * grid.Ca_[mid] *
          ezy_[pml_idx] + com.e_y_coef2_[jt] *
          grid.Cby_[mid] * (grid.hx_[grid.pi(it, jt, kt)] 
                           - grid.hx_[grid.pi(it, jt-1, kt)]);
        
        ezx_[pml_idx] = com.e_x_coef1_[it] * grid.Ca_[mid] *
          ezx_[pml_idx] - com.e_x_coef2_[it] *
          grid.Cbx_[mid] * (grid.hy_[pi(it, jt, kt)] 
                           - grid.hy_[pi(it-1, jt, kt)]);
        
        grid.ez_[grid_idx] = ezx_[pml_idx] + ezy_[pml_idx];
      }
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
