
#include "UPmlCommon.hh"
#include "UPml.hh"
#include "../Types.hh"
#include "../Constants.hh"
#include "../Exceptions.hh"

#include <cmath>

UPmlCommon::UPmlCommon(const Grid &grid)
  : grid_(grid), poly_order_(4), sigma_x_(0), sigma_y_(0), sigma_z_(0), 
    c1_(0), c2_(0), c3_(0), c4_(0), c5_(0), c6_(0), x_size_(0), 
    y_size_(0), z_size_(0)
{
  num_materials_ = grid.get_material_lib()->num_materials();
}

UPmlCommon *UPmlCommon::get_upml_common(Grid &grid)
{
  void *tmp = grid.get_auxdata(UPML_COMMON);

  if (tmp)
    return (UPmlCommon *)tmp;

  UPmlCommon *upml_common = new UPmlCommon(grid);
  upml_common->init_coeffs(grid);

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
  const GridInfo &gi = grid.get_grid_info();

  // Thickesses? 
  for (int i = 0; i < 6; i++)
  {
    thicknesses_[i] = 0;

    if (gi.get_bc_type(static_cast<Face>(i)) == UPML) {
      BoundaryCond &bc = gi.get_boundary(static_cast<Face>(i));
      UPml *p = dynamic_cast<UPml *>(&bc);

      if (p) 
        thicknesses_[i] = p->get_thickness();
    }
  }

  // Set up the sigma array's
  unsigned int szx, szy, szz;
  szx = (thicknesses_[FRONT] + thicknesses_[BACK]) * num_materials_;
  szy = (thicknesses_[LEFT] + thicknesses_[RIGHT]) * num_materials_;
  szz = (thicknesses_[TOP] + thicknesses_[BOTTOM]) * num_materials_;

  sigma_x_ = new float[szx];
  sigma_y_ = new float[szy];
  sigma_z_ = new float[szz];

  if (!sigma_x_ || !sigma_y_ || !sigma_z_)
  {
    free_sigmas();
    throw MemoryException();
  }

  memset(sigma_x_, 0, sizeof(float) * szx);

  memset(sigma_y_, 0, sizeof(float) * szy);

  memset(sigma_z_, 0, sizeof(float) * szz);

  init_sigmas();
}

void UPmlCommon::init_sigmas()
{
  shared_ptr<MaterialLib> mlib = grid_.get_material_lib();
  delta_t delta;
  float *sigmas;
  unsigned int start, thickness;
  int incr;

  map<string, Material>::const_iterator iter = mlib->materials_.begin();
  map<string, Material>::const_iterator iter_e = mlib->materials_.end();
  unsigned int xoffset = 0, yoffset = 0, zoffset = 0;
  
  // Loop over the materials
  while (iter != iter_e) 
  {
    mat_prop_t eps = ((*iter).second).get_epsilon() * EPS_0;
    mat_prop_t sig = ((*iter).second).get_sigma();
    mat_prop_t mu = ((*iter).second).get_mu() * MU_0;
    //mat_prop_t sigs = (*iter).get_sigma_star();

    if (isinf(sig)) 
      //if (sig == __infinity)
    {
      cerr << "UPmlCommon::init_sigmas(): Warning; material is perfect conductor." << endl;
      
    }
    
    // Loop over the faces, for each may have it's own thickness
    for (int faceidx = 0; faceidx < 6; faceidx++)
    {
      switch (faceidx) {
      case FRONT:
      case BACK:
        delta = grid_.get_deltax();
        sigmas = sigma_x_;

        if (faceidx == FRONT)
        {
          start = xoffset + thicknesses_[BACK] + thicknesses_[FRONT] - 1;
          thickness = thicknesses_[FRONT];
          incr = -1;
        } else {
          start = xoffset;
          thickness = thicknesses_[BACK];
          incr = 1;
        }
        break;

      case LEFT:
      case RIGHT:
        delta = grid_.get_deltay();
        sigmas = sigma_y_;

        if (faceidx == RIGHT)
        {
          start = yoffset + thicknesses_[RIGHT] + thicknesses_[LEFT] - 1;
          thickness = thicknesses_[RIGHT];
          incr = -1;
        } else {
          start = yoffset;
          thickness = thicknesses_[LEFT];
          incr = 1;
        }
        break;

      case TOP:
      case BOTTOM:
        delta = grid_.get_deltaz();
        sigmas = sigma_z_;

        if (faceidx == TOP)
        {
          start = zoffset + thicknesses_[TOP] + thicknesses_[BOTTOM] - 1;
          thickness = thicknesses_[TOP];
          incr = -1;
        } else {
          start = zoffset;
          thickness = thicknesses_[BOTTOM];
          incr = 1;
        }
        break;
      }
      mat_coef_t sigma_max = calc_sigma_max(eps, delta);
      
      for (int sigidx = start, idx = thickness; idx > 0; 
           idx--, sigidx += incr)
      {
        sigmas[sigidx] = sigma_max * pow(static_cast<float>(idx) 
					 / static_cast<float>(thickness), 
                                         static_cast<float>(poly_order_));
      }      
    }
    xoffset += thicknesses_[FRONT] + thicknesses_[BACK];
    yoffset += thicknesses_[LEFT] + thicknesses_[RIGHT];
    zoffset += thicknesses_[BOTTOM] + thicknesses_[TOP];
    ++iter;
  }
}
  
mat_coef_t UPmlCommon::calc_sigma_max(mat_prop_t eps, delta_t delta)
{
  return (poly_order_ + 1) / (150 * PI * delta * sqrt(eps));
}
