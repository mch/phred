#include "PmlCommon.hh"
#include "Grid.hh"
#include "Exceptions.hh"

PmlCommon::PmlCommon()
{}

PmlCommon::~PmlCommon()
{}

void PmlCommon::alloc_coeffs(Grid &grid)
{
  free_coeffs();

  ratio_x_ = new float[grid.get_ldx()];
  ratio_star_x_ = new float[grid.get_ldx()];
  
  ratio_y_ = new float[grid.get_ldy()];
  ratio_star_y_ = new float[grid.get_ldy()];
  
  ratio_z_ = new float[grid.get_ldz()];
  ratio_star_z_ = new float[grid.get_ldz()];
  
  e_x_coef1_ = new float[grid.get_ldx()];
  e_x_coef2_ = new float[grid.get_ldx()];
  
  e_y_coef1_ = new float[grid.get_ldy()];
  e_y_coef2_ = new float[grid.get_ldy()];
  
  e_z_coef1_ = new float[grid.get_ldz()];
  e_z_coef2_ = new float[grid.get_ldz()];
  
  h_x_coef1_ = new float[grid.get_ldx()];
  h_x_coef2_ = new float[grid.get_ldx()];
  
  h_y_coef1_ = new float[grid.get_ldy()];
  h_y_coef2_ = new float[grid.get_ldy()];
  
  h_z_coef1_ = new float[grid.get_ldz()];
  h_z_coef2_ = new float[grid.get_ldz()];  

  if (!ratio_x_ || !ratio_star_x_ || !ratio_y_ || !ratio_star_y_ || 
      !ratio_z_ || !ratio_star_z_ || !e_x_coef1_ || !e_x_coef2_ ||
      !e_y_coef1_ || !e_y_coef2_ || !e_z_coef1_ || !e_z_coef2_ ||
      !h_x_coef1_ || !h_x_coef2_ || !h_y_coef1_ || !h_y_coef2_ || 
      !h_z_coef1_ || !h_z_coef2_)
  {
    free_coeffs();
    throw MemoryException();
  }
}

void PmlCommon::free_coeffs()
{
  if (ratio_x_)
    delete[] ratio_x_;

  if (ratio_star_x_)
    delete[] ratio_star_x_;
  
  if (ratio_y_)
    delete[] ratio_y_;

  if (ratio_star_y_)
    delete[] ratio_star_y_;
  
  if (ratio_z_)
    delete[] ratio_z_;
  if (ratio_star_z_)
    delete[] ratio_star_z_;
  
  if (e_x_coef1_)
    delete[] e_x_coef1_;
  if (e_x_coef2_)
    delete[] e_x_coef2_;
  
  if (e_y_coef1_)
    delete[] e_y_coef1_;
  if (e_y_coef2_)
    delete[] e_y_coef2_;
  
  if (e_z_coef1_)
    delete[] e_z_coef1_;
  if (e_z_coef2_)
    delete[] e_z_coef2_;
  
  if (h_x_coef1_)
    delete[] h_x_coef1_;
  if (h_x_coef2_)
    delete[] h_x_coef2_;
  
  if (h_y_coef1_)
    delete[] h_y_coef1_;
  if (h_y_coef2_)
    delete[] h_y_coef2_;
  
  if (h_z_coef1_)
    delete[] h_z_coef1_;
  if (h_z_coef2_)
    delete[] h_z_coef2_;
}


void PmlCommon::init_coeffs(Grid &grid) 
{
  alloc_coeffs(grid);

  GridInfo &gi = grid.get_grid_info();

  for (int i = 0; i < 6; i++)
  {
    if (gi.get_bc_type(static_cast<Face>(i)) == PML) {
      BoundaryCond &bc = gi.get_boundary(static_cast<Face>(i));
      Pml *p = dynamic_cast<Pml *>(&bc);

      if (p) 
      {
        init_ratios(static_cast<Face>(i), grid, p);
      } // if (p)
    }
  } // for
}

void PmlCommon::init_ratios(Face face, Grid &grid, Pml *p)
{
  delta_t d_space;
  int step; 
  unsigned int idx;
  float *arr;

  switch (face)
  {
  case BACK:
    d_space = grid.get_deltax();
    step = 1;
    idx = 0;
    arr  = ratio_x_;
    break;

  case FRONT:
    d_space = grid.get_deltax();
    step = -1; 
    idx = grid.get_ldx() - 1;
    arr  = ratio_x_;
    break;

  case LEFT:
    d_space = grid.get_deltay();
    step = 1;
    idx = 0;
    arr  = ratio_y_;
    break;

  case RIGHT:
    d_space = grid.get_deltay();
    step = -1; 
    idx = grid.get_ldy() - 1;
    arr  = ratio_y_;
    break;

  case BOTTOM:
    d_space = grid.get_deltaz();
    step = 1;
    idx = 0;
    arr  = ratio_z_;
    break;

  case TOP:
    d_space = grid.get_deltaz();
    step = -1; 
    idx = grid.get_ldz() - 1;
    arr  = ratio_z_;
    break;
  }

  for (unsigned int i = 0; i <= p->get_thickness(); i++)
  {
    arr[idx] = 1.0 / d_space 
      * ( p->sigma_over_eps_int((i+0.5)*d_space) 
          - p->sigma_over_eps_int((i-0.5)*d_space));
    
    arr[idx] = 1.0 / d_space 
      * ( p->sigma_over_eps_int((i+1.0)*d_space) 
          - p->sigma_over_eps_int((i)*d_space));
    
    idx += step;
  }

  // And finally, init coefficients
  delta_t dt = grid.get_deltat();

  for(unsigned int i = 0; i < grid.get_ldx(); i++)
  {
    e_x_coef1_[i] = (ratio_x_[i] == 0.0 ? 1.0
                     : exp(-ratio_x_[i] * dt));

    e_x_coef2_[i] = (ratio_x_[i] == 0.0 ? 1.0
                     : (1.0 - exp(-ratio_x_[i] * dt))
                     / (ratio_x_[i] * dt));
                        
    h_x_coef1_[i] = (ratio_star_x_[i] == 0.0 ? 1.0 
                     : exp(-ratio_star_x_[i] * dt));

    h_x_coef2_[i] = (ratio_star_x_[i] == 0.0 ? 1.0 
                     : (1.0 - exp(-ratio_star_x_[i] * dt))
                     / (ratio_star_x_[i] * dt));
  }
                

  for(unsigned int j = 0; j < grid.get_ldy(); j++)
  {
    e_y_coef1_[j] = (ratio_y_[j] == 0.0 ? 1.0 
                     : exp(-ratio_y_[j] * dt));

    e_y_coef2_[j] = (ratio_y_[j] == 0.0 ? 1.0 
                     : (1.0 - exp(-ratio_y_[j] * dt))
                     / (ratio_y_[j] * dt));
                        
    h_y_coef1_[j] = (ratio_star_y_[j] == 0.0 ? 1.0 
                     : exp(-ratio_star_y_[j] * dt));

    h_y_coef2_[j] = (ratio_star_y_[j] == 0.0 ? 1.0 
                     : (1.0 - exp(-ratio_star_y_[j] * dt))
                     / (ratio_star_y_[j] * dt));
  }

  for(unsigned int k = 0; k < grid.get_ldz(); k++)
  {
    e_z_coef1_[k] = (ratio_z_[k] == 0.0 ? 1.0
                     : exp(-ratio_z_[k] * dt));

    e_z_coef2_[k] = (ratio_z_[k] == 0.0 ? 1. 
                     : (1.0 - exp(-ratio_z_[k] * dt))
                     / (ratio_z_[k] * dt));
                        
    h_z_coef1_[k] = (ratio_star_z_[k] == 0.0 ? 1.0 
                     : exp(-ratio_star_z_[k] * dt));

    h_z_coef2_[k] = (ratio_star_z_[k] == 0.0 ? 1.0 
                     : (1.0 - exp(-ratio_star_z_[k] * dt))
                     / (ratio_star_z_[k] * dt));
  }
                  

}

