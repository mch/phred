
#include "UPmlCommon.hh"

UPmlCommon *UPmlCommon::get_upml_common(Grid &grid)
{
  if (grid.get_define_mode())
    return 0;

  void *tmp = grid.get_auxdata(UPML_COMMON);

  if (tmp)
    return (UpmlCommon *)tmp;

  UPmlCommon *upml_common = new UPmlCommon(grid);
  upml_common->init_coeffs();

  grid.add_auxdata(UPML_COMMON, (void *)upml_common);

  return upml_common;
}

void UPmlCommon::free_sigmas()
{
  if (sigma_x_)
  {
    delete[] sigma_x_;
    sigma_x_ = 0;
  }

  if (sigma_y_)
  {
    delete[] sigma_y_;
    sigma_y_ = 0;
  }

  if (sigma_z_)
  {
    delete[] sigma_z_;
    sigma_z_ = 0;
 }
}

void UPmlCommon::init_coeffs(Grid &grid)
{
  GridInfo &gi = grid.get_grid_info();

  // Set up the sigma array's
  sigma_x_ = new float[grid.get_ldx()];
  sigma_y_ = new float[grid.get_ldy()];
  sigma_z_ = new float[grid.get_ldz()];

  if (!sigma_x_ || !sigma_y_ || !sigma_z_)
  {
    free_sigmas();
    throw MemoryException();
  }

  memset(sigma_x_, 0, sizeof(float) * grid.get_ldx());
  memset(sigma_y_, 0, sizeof(float) * grid.get_ldy());
  memset(sigma_z_, 0, sizeof(float) * grid.get_ldz());

  for (int i = 0; i < 6; i++)
  {
    if (gi.get_bc_type(static_cast<Face>(i)) == PML) {
      BoundaryCond &bc = gi.get_boundary(static_cast<Face>(i));
      UPml *p = dynamic_cast<UPml *>(&bc);

      if (p) 
      {
        init_sigmas(static_cast<Face>(i), grid, p);
      } // if (p)
    }
  } // for

}

void UPmlCommon::init_sigmas(Face face, const Grid &grid, UPml *p)
{
  delta_t d_space;
  int step; 
  unsigned int idx1, idx2;
  float *arr1, *arr2;

  switch (face)
  {
  case BACK:
    d_space = grid.get_deltax();
    step = 1;
    idx1 = 0;
    idx2 = 0;
    arr1 = ratio_x_;
    arr2 = ratio_star_x_;
    break;

  case FRONT:
    d_space = grid.get_deltax();
    step = -1; 
    idx1 = grid.get_ldx_sd() - 1;
    idx2 = grid.get_ldx_sd() - 2;
    arr1 = ratio_x_;
    arr2 = ratio_star_x_;
    break;

  case LEFT:
    d_space = grid.get_deltay();
    step = 1;
    idx1 = 0;
    idx2 = 0;
    arr1 = ratio_y_;
    arr2 = ratio_star_y_;
    break;

  case RIGHT:
    d_space = grid.get_deltay();
    step = -1; 
    idx1 = grid.get_ldy_sd() - 1;
    idx2 = grid.get_ldy_sd() - 2;
    arr1 = ratio_y_;
    arr2 = ratio_star_y_;
    break;

  case BOTTOM:
    d_space = grid.get_deltaz();
    step = 1;
    idx1 = 0;
    idx2 = 0;
    arr1 = ratio_z_;
    arr2 = ratio_star_z_;
    break;

  case TOP:
    d_space = grid.get_deltaz();
    step = -1; 
    idx1 = grid.get_ldz_sd() - 1;
    idx2 = grid.get_ldz_sd() - 2;
    arr1 = ratio_z_;
    arr2 = ratio_star_z_;
    break;
  }

  mat_coef_t sigma_max = (poly_order_ + 1) 
    / (150 * PI * d_space * sqrt((*iter).get_epsilon()));

  for (unsigned int i = 0; i <= p->get_thickness(); i++)
  {
    arr1[idx1] = 1.0 / d_space 
      * ( p->sigma_over_eps_int((i+0.5)*d_space) 
          - p->sigma_over_eps_int((i-0.5)*d_space));
    
    arr2[idx2] = 1.0 / d_space 
      * ( p->sigma_over_eps_int((i+1.0)*d_space) 
          - p->sigma_over_eps_int((i)*d_space));
    
    idx1 += step;
    idx2 += step;
  }

}
