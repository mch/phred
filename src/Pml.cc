#include "Pml.hh"
#include "Grid.hh"
#include "Exceptions.hh"
#include "Constants.hh"

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

  unsigned int sz = (r.xmax - r.xmin) * (r.ymax - r.ymin) 
    * (r.zmax - r.zmin);
  
  sz_ = sz;

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
    ratio_m_ *= log(g_) / (pow(g_, static_cast<float>(thickness_)) - 1.0);
    geometric_profile_ = 1;
    break;
  }

//   cout << "Pml setup TEST results for face " << face << "\n" 
//        << "\tratio_m_ = " << ratio_m_
//        << "\n\texponent_n_ = " << exponent_n_
//        << "\n\tdelta_bndy_ = " << delta_bndy_
//        << "\n\tgeometric_delta_ = " << geometric_delta_
//        << "\n\tgeomtric_profile_ = " << geometric_profile_ << endl;

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

//   if (eyx_[12055])
//     cout << "1) EYX_ HAS BEEN CORRUPTED!!!" << endl;
//   else 
//     cout << "1) eyx_ ok." << endl;

  if (type == E) {

    region_t e_grid_r = grid_r; 
    region_t e_pml_r = pml_r_;
    
    // Modify the grid region so that the outer walls are not
    // computed; that they be electric walls.
    if (grid_r.xmin == 0) 
    {
      e_grid_r.xmin++;
      e_pml_r.xmin++;
    }

    if (grid_r.ymin == 0)
    {
      e_grid_r.ymin++;
      e_pml_r.ymin++;
    }

    if (grid_r.zmin == 0)
    {
      e_grid_r.zmin++;
      e_pml_r.zmin++;
    }

    if (grid_r.xmax == grid.get_ldx())
    {
      e_grid_r.xmax--;
      e_pml_r.xmax--;
    }

    if (grid_r.ymax == grid.get_ldy())
    {
      e_grid_r.ymax--;
      e_pml_r.ymax--;
    }

    if (grid_r.zmax == grid.get_ldz())
    {
      e_grid_r.zmax--;
      e_pml_r.zmax--;
    }

//      cout << "Electric field Pml update on face " << static_cast<int>(face) 
//           << " in grid ranges x={" << e_grid_r.xmin << "," 
//          << e_grid_r.xmax << "}, y={" << e_grid_r.ymin << "," 
//          << e_grid_r.ymax << "}, z={" << e_grid_r.zmin << ","
//          << e_grid_r.zmax << "}.\n"
//          << "pml range x={" << e_pml_r.xmin << "," 
//          << e_pml_r.xmax << "}, y={" << e_pml_r.ymin << "," 
//          << e_pml_r.ymax << "}, z={" << e_pml_r.zmin << ","
//          << e_pml_r.zmax << "}.\n";

    pml_update_ex(e_pml_r, e_grid_r, grid_r, grid);
//     if (eyx_[12055])
//       cout << "2) EYX_ HAS BEEN CORRUPTED!!!" << endl;
//     else 
//       cout << "2) eyx_ ok." << endl;
    
    pml_update_ey(e_pml_r, e_grid_r, grid_r, grid);
//     if (eyx_[12055])
//       cout << "3) EYX_ HAS BEEN CORRUPTED!!!" << endl;
//     else 
//       cout << "3) eyx_ ok." << endl;

    pml_update_ez(e_pml_r, e_grid_r, grid_r, grid);
//     if (eyx_[12055])
//       cout << "4) EYX_ HAS BEEN CORRUPTED!!!" << endl;
//     else 
//       cout << "4) eyx_ ok." << endl;

  }
  else if (type == H)
  {
//     cout << "Magnetic field Pml update on face " << static_cast<int>(face) 
//          << " in grid ranges x={" << grid_r.xmin << "," 
//          << grid_r.xmax << "}, y={" << grid_r.ymin << "," 
//          << grid_r.ymax << "}, z={" << grid_r.zmin << ","
//          << grid_r.zmax << "}.\n"
//          << "pml range x={" << pml_r_.xmin << "," 
//          << pml_r_.xmax << "}, y={" << pml_r_.ymin << "," 
//          << pml_r_.ymax << "}, z={" << pml_r_.zmin << ","
//          << pml_r_.zmax << "}.\n";

    pml_update_hx(grid_r, grid);
    pml_update_hy(grid_r, grid);
    pml_update_hz(grid_r, grid);
  } 
  else
  {
    cout << "INCORRECT FIELD TYPE GIVEN TO UPDATE!" << endl;
  }
  
}

void Pml::pml_update_ex(const region_t &e_pml_r, 
                        const region_t &e_grid_r, 
                        const region_t &grid_r, 
                        Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  PmlCommon &com = grid.get_pml_common();

  for(i = 0, it = grid_r.xmin; it < grid_r.xmax - 1; i++, it++)
    for(j = e_pml_r.ymin, jt = e_grid_r.ymin; jt < e_grid_r.ymax; j++, jt++)
      for(k = e_pml_r.zmin, kt = e_grid_r.zmin; kt < e_grid_r.zmax; k++, kt++)
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
        
//         if (grid.hy_[grid.pi(it, jt, kt)] 
//             - grid.hy_[grid.pi(it, jt, kt-1)] > 0.0  || 
//             grid.hz_[grid.pi(it, jt, kt)] 
//             - grid.hz_[grid.pi(it, jt-1, kt)] > 0.0)
//         {
//           cout << "ex pml calc at " << it << "," << jt << "," << kt 
//                << "; exz_ = " << exz_[pml_idx] << ", exy_ = "
//                << exy_[pml_idx] << endl;
//         }

        grid.ex_[grid_idx] = exz_[pml_idx] + exy_[pml_idx];
      }
}

void Pml::pml_update_ey(const region_t &e_pml_r, 
                        const region_t &e_grid_r, 
                        const region_t &grid_r, 
                        Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  PmlCommon &com = grid.get_pml_common();

  for(i = e_pml_r.xmin, it = e_grid_r.xmin; it < e_grid_r.xmax; i++, it++)
    for(j = 0, jt = grid_r.ymin; jt < grid_r.ymax - 1; j++, jt++)
      for(k = e_pml_r.zmin, kt = e_grid_r.zmin; kt < e_grid_r.zmax; k++, kt++)
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
        
//         if (grid.hz_[grid.pi(it, jt, kt)] 
//             - grid.hz_[grid.pi(it-1, jt, kt)] > 0.0  || 
//             grid.hx_[grid.pi(it, jt, kt)] 
//             - grid.hx_[grid.pi(it, jt, kt-1)] > 0.0)
//         {
//           cout << "ey pml calc at " << it << "," << jt << "," << kt 
//                << "; eyx_ = " << eyx_[pml_idx] << ", eyz_ = " 
//                << eyz_[pml_idx] << endl;
//         }

        grid.ey_[grid_idx] = eyx_[pml_idx] + eyz_[pml_idx];
      }
}

void Pml::pml_update_ez(const region_t &e_pml_r, 
                        const region_t &e_grid_r, 
                        const region_t &grid_r, 
                        Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  PmlCommon &com = grid.get_pml_common();

  for(i = e_pml_r.xmin, it = e_grid_r.xmin; it < e_grid_r.xmax; i++, it++)
    for(j = e_pml_r.ymin, jt = e_grid_r.ymin; jt < e_grid_r.ymax; j++, jt++)
      for(k = 0, kt = grid_r.zmin; kt < grid_r.zmax - 1; k++, kt++)
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
        
//         if (grid.hy_[grid.pi(it, jt, kt)] 
//             - grid.hy_[grid.pi(it-1, jt, kt)] > 0.0  || 
//             grid.hx_[grid.pi(it, jt, kt)] 
//             - grid.hx_[grid.pi(it, jt-1, kt)] > 0.0)
//         {
//           cout << "ez pml calc at " << it << "," << jt << "," << kt 
//                << "; ezy_ = " << ezy_[pml_idx] << ", ezx_ = "
//                << ezx_[pml_idx] << endl;
//         }

        grid.ez_[grid_idx] = ezx_[pml_idx] + ezy_[pml_idx];
      }
}

void Pml::pml_update_hx(const region_t &grid_r, Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  PmlCommon &com = grid.get_pml_common();

  for(i = pml_r_.xmin, it = grid_r.xmin; it < grid_r.xmax; i++, it++)
    for(j = pml_r_.ymin, jt = grid_r.ymin; jt < grid_r.ymax-1; j++, jt++)
      for(k = pml_r_.zmin, kt = grid_r.zmin; kt < grid_r.zmax-1; k++, kt++)
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
        
//         if (grid.ey_[grid.pi(it, jt, kt+1)] 
//             - grid.ey_[grid.pi(it, jt, kt)] > 0.0  || 
//             grid.ez_[grid.pi(it, jt+1, kt)] 
//             - grid.ez_[grid.pi(it, jt, kt)] > 0.0)
//         {
//           cout << "hx pml calc at " << it << "," << jt << "," << kt 
//                << "; hxz_ = " << hxz_[pml_idx] << ", hxy_ = " 
//                << hxy_[pml_idx] << endl;
//         }

        grid.hx_[grid_idx] = hxz_[pml_idx] + hxy_[pml_idx];
      }
}

void Pml::pml_update_hy(const region_t &grid_r, Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  PmlCommon &com = grid.get_pml_common();

  for(i = pml_r_.xmin, it = grid_r.xmin; it < grid_r.xmax - 1; i++, it++)
    for(j = pml_r_.ymin, jt = grid_r.ymin; jt < grid_r.ymax; j++, jt++)
      for(k = pml_r_.zmin, kt = grid_r.zmin; kt < grid_r.zmax-1; k++, kt++)
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
        
//         if (grid.ez_[grid.pi(it+1, jt, kt)] 
//             - grid.ez_[grid.pi(it, jt, kt)] > 0.0  || 
//             grid.ex_[grid.pi(it, jt, kt+1)] 
//             - grid.ex_[grid.pi(it, jt, kt)] > 0.0)
//         {
//           cout << "hy pml calc at " << it << "," << jt << "," << kt 
//                << "; hyx_ = " << hyx_[pml_idx] << ", hyz_ = "
//                << hyz_[pml_idx] << endl;
//         }

        grid.hy_[grid_idx] = hyx_[pml_idx] + hyz_[pml_idx];
      }
}

void Pml::pml_update_hz(const region_t &grid_r, Grid &grid)
{
  unsigned int grid_idx, pml_idx, mid; 

  unsigned int i,j,k; 	/* indices in PML-layer */
  unsigned int it,jt,kt;/* indices in total computational domain (FDTD grid) */

  PmlCommon &com = grid.get_pml_common();

  for(i = pml_r_.xmin, it = grid_r.xmin; it < grid_r.xmax-1; i++, it++)
    for(j = pml_r_.ymin, jt = grid_r.ymin; jt < grid_r.ymax-1; j++, jt++)
      for(k = pml_r_.zmin, kt = grid_r.zmin; kt < grid_r.zmax; k++, kt++)
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
        
//         if (grid.ex_[grid.pi(it, jt+1, kt)] 
//             - grid.ex_[grid.pi(it, jt, kt)] > 0.0  || 
//             grid.ey_[grid.pi(it+1, jt, kt)] 
//             - grid.ey_[grid.pi(it, jt, kt)] > 0.0)
//         {
//           cout << "hz pml calc at " << it << "," << jt << "," << kt 
//                << "; hzy_ = " << hzy_[pml_idx] << ", hzx_ = "
//                << hzx_[pml_idx] << endl;
//         }

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

RxTxData Pml::get_rx_tx_data(Face pmlface, Face sdface)
{
  RxTxData ret;

  switch (sdface) {
  case BACK:
    break;
  case FRONT:
    break;
  case LEFT:
    break;
  case RIGHT:
    break;
  case BOTTOM:
    break;
  case TOP:
    break;
  }

  return ret;
}
