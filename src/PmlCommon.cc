#include "PmlCommon.hh"
#include "Grid.hh"

PmlCommon::PmlCommon()
{}

PmlCommon::~PmlCommon()
{}

void PmlCommon::alloc_coeffs(Grid &grid)
{
  free_coeffs();

  ratio_x_ = new float[grid.get_gdx()];
  ratio_star_x_ = new float[grid.get_gdx()];
  
  ratio_y_ = new float[grid.get_gdy()];
  ratio_star_y_ = new float[grid.get_gdy()];
  
  ratio_z_ = new float[grid.get_gdz()];
  ratio_star_z_ = new float[grid.get_gdz()];
  
  e_x_coef1_ = new float[grid.get_gdx()];
  e_x_coef2_ = new float[grid.get_gdx()];
  
  e_y_coef1_ = new float[grid.get_gdy()];
  e_y_coef2_ = new float[grid.get_gdy()];
  
  e_z_coef1_ = new float[grid.get_gdz()];
  e_z_coef2_ = new float[grid.get_gdz()];
  
  h_x_coef1_ = new float[grid.get_gdx()];
  h_x_coef2_ = new float[grid.get_gdx()];
  
  h_y_coef1_ = new float[grid.get_gdy()];
  h_y_coef2_ = new float[grid.get_gdy()];
  
  h_z_coef1_ = new float[grid.get_gdz()];
  h_z_coef2_ = new float[grid.get_gdz()];  
}

void PmlCommon::free_coeffs()
{
  if (!ratio_x_)
    return;
  
  delete[] ratio_x_;
  delete[] ratio_x_;
  delete[] ratio_star_x_;
  
  delete[] ratio_y_;
  delete[] ratio_star_y_;
  
  delete[] ratio_z_;
  delete[] ratio_star_z_;
  
  delete[] e_x_coef1_;
  delete[] e_x_coef2_;
  
  delete[] e_y_coef1_;
  delete[] e_y_coef2_;
  
  delete[] e_z_coef1_;
  delete[] e_z_coef2_;
  
  delete[] h_x_coef1_;
  delete[] h_x_coef2_;
  
  delete[] h_y_coef1_;
  delete[] h_y_coef2_;
  
  delete[] h_z_coef1_;
  delete[] h_z_coef2_;
}


void PmlCommon::init_coeffs(Grid &grid) 
{
  GridInfo &gi = grid.get_grid_info();

  for (int i = 0; i < 6; i++)
  {
    if (gi.get_bc_type(static_cast<Face>(i)) == PML) {
      BoundaryCond &bc = gi.get_boundary(static_cast<Face>(i));
      
      // initialize maximum ratios ...
    }
  }
}
