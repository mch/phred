
#include "UPmlCommon.hh"
#include "UPml.hh"
#include "../Types.hh"
#include "../Constants.hh"
#include "../Exceptions.hh"

#include <cmath>

UPmlCommon::UPmlCommon(const Grid &grid)
  : grid_(grid), inited_(false), poly_order_(4), 
    sigma_max_(0.0), eps_opt_(1.0), 
    sigma_ratio_(1.0),
    x_size_(0), y_size_(0), z_size_(0),
    sigma_x_(0), sigma_y_(0), sigma_z_(0), 
    kx_(0), ky_(0), kz_(0),
    Ax_(0), Ay_(0), Az_(0), 
    Bx_(0), By_(0), Bz_(0), 
    Cx_(0), Cy_(0), Cz_(0), 
    Dx_(0), Dy_(0), Dz_(0), 
    er_(0), ur_(0)
{

}

UPmlCommon *UPmlCommon::get_upml_common(Grid &grid)
{
  void *tmp = grid.get_auxdata(UPML_COMMON);

  if (tmp)
    return (UPmlCommon *)tmp;

  UPmlCommon *upml_common = new UPmlCommon(grid);

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

  if (kx_)
  {
    delete[] kx_;
    kx_ = 0;
    delete[] ky_;
    ky_ = 0;
    delete[] kz_;
    kz_ = 0;
  }
}

void UPmlCommon::init_coeffs()
{
  if (inited_)
    return;
  
  const GridInfo &gi = grid_.get_grid_info();

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
//   szx = thicknesses_[FRONT] + thicknesses_[BACK];
//   szy = thicknesses_[LEFT] + thicknesses_[RIGHT];
//   szz = thicknesses_[TOP] + thicknesses_[BOTTOM];

  x_size_ = szx = grid_.get_ldx_sd();
  y_size_ = szy = grid_.get_ldy_sd();
  z_size_ = szz = grid_.get_ldz_sd();

  // Conductivity
  sigma_x_ = new float[szx];
  sigma_y_ = new float[szy];
  sigma_z_ = new float[szz];

  kx_ = new float[x_size_];
  ky_ = new float[y_size_];
  kz_ = new float[z_size_];

  Ax_ = new float[x_size_];
  Ay_ = new float[y_size_];
  Az_ = new float[z_size_];

  Bx_ = new float[x_size_];
  By_ = new float[y_size_];
  Bz_ = new float[z_size_];

  Cx_ = new float[x_size_];
  Cy_ = new float[y_size_];
  Cz_ = new float[z_size_];

  Dx_ = new float[x_size_];
  Dy_ = new float[y_size_];
  Dz_ = new float[z_size_];

  memset(sigma_x_, 0, sizeof(float) * szx);
  memset(sigma_y_, 0, sizeof(float) * szy);
  memset(sigma_z_, 0, sizeof(float) * szz);

  memset(kx_, 0, sizeof(float) * x_size_);
  memset(ky_, 0, sizeof(float) * y_size_);
  memset(kz_, 0, sizeof(float) * z_size_);

  for (int i = 0; i < x_size_; i++)
    kx_[i] = 1;

  for (int i = 0; i < y_size_; i++)
    ky_[i] = 1;

  for (int i = 0; i < z_size_; i++)
    kz_[i] = 1;

  memset(Ax_, 0, sizeof(float) * x_size_);
  memset(Ay_, 0, sizeof(float) * y_size_);
  memset(Az_, 0, sizeof(float) * z_size_);

  memset(Bx_, 0, sizeof(float) * x_size_);
  memset(By_, 0, sizeof(float) * y_size_);
  memset(Bz_, 0, sizeof(float) * z_size_);

  memset(Cx_, 0, sizeof(float) * x_size_);
  memset(Cy_, 0, sizeof(float) * y_size_);
  memset(Cz_, 0, sizeof(float) * z_size_);

  memset(Dx_, 0, sizeof(float) * x_size_);
  memset(Dy_, 0, sizeof(float) * y_size_);
  memset(Dz_, 0, sizeof(float) * z_size_);

  int nm = grid_.get_material_lib()->num_materials();

  er_ = new float[nm];
  ur_ = new float[nm];
  mtype_ = new MaterialType[nm];
  lossyA_ = new float[nm];
  lossyB_ = new float[nm];

  memset(er_, 0, sizeof(float) * nm);
  memset(ur_, 0, sizeof(float) * nm);
  memset(mtype_, 0, sizeof(MaterialType) * nm);
  memset(lossyA_, 0, sizeof(float) * nm);
  memset(lossyB_, 0, sizeof(float) * nm);

  init_sigmas();

  init_constants();

  inited_ = true;
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
  
  // Loop over the materials
  while (iter != iter_e) 
  {
    mat_prop_t eps = ((*iter).second).get_epsilon() * EPS_0;
    mat_prop_t sig = ((*iter).second).get_sigma();
    mat_prop_t mu = ((*iter).second).get_mu() * MU_0;
    //mat_prop_t sigs = (*iter).get_sigma_star();

    mat_idx_t mid = (*iter).second.get_id();

    er_[mid] = 1 / (eps);
    ur_[mid] = 1 / (mu);
    mtype_[mid] = (*iter).second.type();

    if (mtype_[mid] == LOSSY)
    {
      delta_t dt = grid_.get_deltat();
      lossyA_[mid] = (2*eps - dt*sig) / (2*eps + dt*sig);
      lossyB_[mid] = (2*dt) / (2*eps + dt*sig);
    }

    ++iter;
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
        start = grid_.get_ldx_sd() - 2; 
        thickness = thicknesses_[FRONT];
        incr = -1;
      } else {
        start = 0;
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
        start = grid_.get_ldy_sd() - 2;
        thickness = thicknesses_[RIGHT];
        incr = -1;
      } else {
        start = 0;
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
        start = grid_.get_ldz_sd() - 2;
        thickness = thicknesses_[TOP];
        incr = -1;
      } else {
        start = 0;
        thickness = thicknesses_[BOTTOM];
        incr = 1;
      }
      break;
    }

    mat_coef_t sigma_max = calc_sigma_max(delta);
    //mat_coef_t k_max = calc_k_max(delta);

    for (int sigidx = start, idx = thickness; idx > 0; 
         idx--, sigidx += incr)
    {
      sigmas[sigidx] = sigma_max * pow(static_cast<float>(idx) 
                                       / static_cast<float>(thickness), 
                                       static_cast<float>(poly_order_));
    }      
  }

#ifdef DEBUG
  cout << "UPML conductivity profiles: \nx axis: ";

  for (int i = 0; i < grid_.get_ldx_sd(); i++)
    cout << sigma_x_[i] << " ";

  cout << "\ny axis: ";

  for (int i = 0; i < grid_.get_ldy_sd(); i++)
    cout << sigma_y_[i] << " ";

  cout << "\nz axis: ";

  for (int i = 0; i < grid_.get_ldz_sd(); i++)
    cout << sigma_z_[i] << " ";
  
  cout << endl;
#endif

}

void UPmlCommon::init_constants()
{
  for (int i = 0; i < grid_.get_ldx_sd(); i++)
  {
    Ax_[i] = ( 2 * EPS_0 * kx_[i] - grid_.get_deltat() * sigma_x_[i] ) 
      / ( 2 * EPS_0 * kx_[i] + grid_.get_deltat() * sigma_x_[i] );

    Bx_[i] = 2 * EPS_0 * grid_.get_deltat()
      / ( 2 * EPS_0 * kx_[i] + grid_.get_deltat() * sigma_x_[i] )
      ;

    Cx_[i] = ( 2 * EPS_0 * kx_[i] + grid_.get_deltat() * sigma_x_[i] ) 
      / (2 * EPS_0 * grid_.get_deltat());

    Dx_[i] = ( 2 * EPS_0 * kx_[i] - grid_.get_deltat() * sigma_x_[i] ) 
      / (2 * EPS_0 * grid_.get_deltat());
  }

  for (int j = 0; j < grid_.get_ldy_sd(); j++)
  {
    Ay_[j] = ( 2 * EPS_0 * ky_[j] - grid_.get_deltat() * sigma_y_[j] ) 
      / ( 2 * EPS_0 * ky_[j] + grid_.get_deltat() * sigma_y_[j] );

    By_[j] = 2 * EPS_0 * grid_.get_deltat()
      / ( 2 * EPS_0 * ky_[j] + grid_.get_deltat() * sigma_y_[j] )
      ;

    Cy_[j] = ( 2 * EPS_0 * ky_[j] + grid_.get_deltat() * sigma_y_[j] ) 
      / (2 * EPS_0 * grid_.get_deltat());

    Dy_[j] = ( 2 * EPS_0 * ky_[j] - grid_.get_deltat() * sigma_y_[j] ) 
      / (2 * EPS_0 * grid_.get_deltat());
  }

  for (int k = 0; k < grid_.get_ldz_sd(); k++)
  {
    Az_[k] = ( 2 * EPS_0 * kz_[k] - grid_.get_deltat() * sigma_z_[k] ) 
      / ( 2 * EPS_0 * kz_[k] + grid_.get_deltat() * sigma_z_[k] );

    Bz_[k] = 2 * EPS_0 * grid_.get_deltat()
      / ( 2 * EPS_0 * kz_[k] + grid_.get_deltat() * sigma_z_[k] )
      ;

    Cz_[k] = ( 2 * EPS_0 * kz_[k] + grid_.get_deltat() * sigma_z_[k] ) 
      / (2 * EPS_0 * grid_.get_deltat());

    Dz_[k] = ( 2 * EPS_0 * kz_[k] - grid_.get_deltat() * sigma_z_[k] ) 
      / (2 * EPS_0 * grid_.get_deltat());

  }

}

mat_coef_t UPmlCommon::calc_sigma_max(delta_t delta)
{
  mat_coef_t ret = 0.0;
  // Bad floating point test!
  if (sigma_max_ == 0.0)
  {
    ret = sigma_ratio_ * ((poly_order_ + 1) 
                          / (150 * PI * delta * sqrt(eps_opt_)));
#ifdef DEBUG
    cout << "Calculated Maximum conductivity using poly_order = " 
         << poly_order_ << ", eps_opt = " << eps_opt_ 
         << "\nand sigma_ratio = " << sigma_ratio_ << ": " << ret << endl;
#endif
  }
  else
    ret = sigma_max_;

  return ret;
}
