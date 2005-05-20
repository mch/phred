
#include "UPmlCommon.hh"
#include "UPml.hh"
#include "../Types.hh"
#include "../Constants.hh"
#include "../Exceptions.hh"

#include <cmath>

UPmlCommon::UPmlCommon(const Grid &grid)
  : grid_(grid), inited_(false), poly_order_(4), 
    sigma_max_(0.0), eps_opt_(1.0), 
    sigma_ratio_(1.0), k_max_(1),
    x_size_(0), y_size_(0), z_size_(0),
    sigma_x_(0), sigma_y_(0), sigma_z_(0), 
    kx_(0), ky_(0), kz_(0),
    Ax_(0), Ay_(0), Az_(0), 
    Bx_(0), By_(0), Bz_(0), 
    Cx_(0), Cy_(0), Cz_(0), 
    Dx_(0), Dy_(0), Dz_(0), 
    er_(0), ur_(0), mtype_(0),
    lossyA_(0), lossyB_(0),
#ifdef ADE_DRUDE
    drudeC1_(0), drudeC2_(0), drudeC3_(0), drudeC4_(0), drudeC5_(0),
#else
    vcdt_(0), omegasq_(0), eps_inf_(0),
#endif
    debyeA_(0), debyeB_(0), debyeC_(0)
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

  memset(er_, 0, sizeof(float) * nm);
  memset(ur_, 0, sizeof(float) * nm);
  memset(mtype_, 0, sizeof(MaterialType) * nm);

  init_sigmas();

  init_constants();

  inited_ = true;
}

void UPmlCommon::init_sigmas()
{
  shared_ptr<MaterialLib> mlib = grid_.get_material_lib();
  delta_t delta;
  float *sigmas;
  float *k;
  unsigned int start, thickness;
  int incr;
  int nm = grid_.get_material_lib()->num_materials();

  map<string, Material>::const_iterator iter = mlib->materials_.begin();
  map<string, Material>::const_iterator iter_e = mlib->materials_.end();
  
  // Check if despersion constants and aux variables should be allocated
  bool need_lossy = false;
  bool need_drude = false;
  bool need_debye = false;

  while (iter != iter_e) 
  {
    mat_idx_t mid = (*iter).second.get_id();
  
    if ((*iter).second.type() == LOSSY)
      need_lossy = true;

    if ((*iter).second.type() == DRUDE)
      need_drude = true;
    
    if ((*iter).second.type() == DEBYE)
      need_debye = true;

    ++iter;
  }

  if (need_lossy)
  {
    lossyA_ = new float[nm];
    lossyB_ = new float[nm];

    memset(lossyA_, 0, sizeof(float) * nm);
    memset(lossyB_, 0, sizeof(float) * nm);
  }

  if (need_drude)
  {
#ifdef ADE_DRUDE
    drudeC1_ = new float[nm];
    drudeC2_ = new float[nm];
    drudeC3_ = new float[nm];
    drudeC4_ = new float[nm];
    drudeC5_ = new float[nm];

    memset(drudeC1_, 0, sizeof(float) * nm);
    memset(drudeC2_, 0, sizeof(float) * nm);
    memset(drudeC3_, 0, sizeof(float) * nm);
    memset(drudeC4_, 0, sizeof(float) * nm);
    memset(drudeC5_, 0, sizeof(float) * nm);
#else
    vcdt_ = new mat_prop_t[nm];
    omegasq_ = new mat_prop_t[nm];
    eps_inf_ = new mat_prop_t[nm];
#endif
  }

  if (need_debye)
  {
    debyeA_ = new field_t[nm];
    debyeB_ = new field_t[nm];
    debyeC_ = new field_t[nm];
    memset(debyeA_, 0, sizeof(field_t) * nm);
    memset(debyeB_, 0, sizeof(field_t) * nm);
    memset(debyeC_, 0, sizeof(field_t) * nm);
  }

  // Loop over the materials
  iter = mlib->materials_.begin();
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

    if (mtype_[mid] == DRUDE)
    {
      delta_t dt = grid_.get_deltat();
      mat_prop_t omega_p = (*iter).second.get_plasma_freq();
      mat_prop_t vc = (*iter).second.get_collision_freq();

      try {
        eps_inf_[mid] = (iter->second).get_property("drude_epsilon_inf");
      } 
      catch (MaterialPropertyException e)
      {
        eps_inf_[mid] = 0.0;
      }

#ifdef ADE_DRUDE
      mat_prop_t drudeC6 = eps * (
                                  (omega_p * omega_p / 4) 
                                  + (vc / (2 * dt))
                                  - (1 / (dt * dt))
                                  );

      drudeC1_[mid] = eps / drudeC6 
        * ( -(omega_p * omega_p) / 2 - 2 / (dt * dt));
                                             
      drudeC2_[mid] = eps / drudeC6
        * ( 1 / (dt * dt) + vc / (2 * dt) - (omega_p * omega_p) / 4 );
        
      drudeC3_[mid] = 1 / drudeC6 
        * ( vc / (2 * dt) - 1 / (dt * dt));

      drudeC4_[mid] = 2 / (drudeC6 * dt * dt);

      drudeC5_[mid] = (-vc * dt - 2) 
        / (2 * dt * dt * drudeC6);
#else
// #ifdef ISHIMARU_DRUDE
      vcdt_[mid] = exp(-1.0 * vc * dt);
// #else
//       vcdt_[mid] = exp(vc * dt);
// #endif
      omegasq_[mid] = omega_p * omega_p * (dt / vc);
#endif
    }

    if (mtype_[mid] == DEBYE) 
    {
      delta_t dt = grid_.get_deltat();
      mat_prop_t eps_inf = (*iter).second.get_property("debye_eps_inf");
      mat_prop_t eps_s = (*iter).second.get_property("debye_eps_s");
      mat_prop_t tau = (*iter).second.get_property("debye_tau");
      
      debyeA_[mid] = (2 * eps_inf * tau 
                        - eps_s * dt)
        / (2 * eps_inf * tau 
           + eps_s * dt);

      debyeB_[mid] = tau / dt + 0.5;
      debyeC_[mid] = tau / dt - 0.5;      
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
      k = kx_;

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
      k = ky_;

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
      k = kz_;

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

      // Assign k... but how? 
      k[sigidx] = 1 + (k_max_ - 1) * pow(static_cast<float>(idx) 
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

  cout << "UPML k profiles: \nx axis: ";

  for (int i = 0; i < grid_.get_ldx_sd(); i++)
    cout << kx_[i] << " ";

  cout << "\ny axis: ";

  for (int i = 0; i < grid_.get_ldy_sd(); i++)
    cout << ky_[i] << " ";

  cout << "\nz axis: ";

  for (int i = 0; i < grid_.get_ldz_sd(); i++)
    cout << kz_[i] << " ";
  
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
    cout << "UPML: Calculated Maximum conductivity using poly_order = " 
         << poly_order_ << ", eps_opt = " << eps_opt_ 
         << "\nand sigma_ratio = " << sigma_ratio_ << ": " << ret << endl;
#endif
  }
  else
    ret = sigma_max_;

  return ret;
}


void UPmlCommon::deinit()
{
  free_sigmas();

  if (Ax_)
  {
    delete[] Ax_;
    Ax_ = 0;
    delete[] Ay_;
    Ay_ = 0;
    delete[] Az_;
    Az_ = 0;

    delete[] Bx_;
    Bx_ = 0;
    delete[] By_;
    By_ = 0;
    delete[] Bz_;
    Bz_ = 0;

    delete[] Cx_;
    Cx_ = 0;
    delete[] Cy_;
    Cy_ = 0;
    delete[] Cz_;
    Cz_ = 0;

    delete[] Dx_;
    Dx_ = 0;
    delete[] Dy_;
    Dy_ = 0;
    delete[] Dz_;
    Dz_ = 0;
  }

  if (er_)
  {
    delete[] er_;
    er_ = 0;

    delete[] ur_;
    ur_ = 0;

    delete[] mtype_;
    mtype_ = 0;
  }

  if (lossyA_)
  {
    delete[] lossyA_;
    lossyA_ = 0;

    delete[] lossyB_;
    lossyB_ = 0;
  }

#ifdef ADE_DRUDE
  if (drudeC1_)
  {
    delete[] drudeC1_;
    drudeC1_ = 0;

    delete[] drudeC2_;
    drudeC2_ = 0;

    delete[] drudeC3_;
    drudeC3_ = 0;

    delete[] drudeC4_;
    drudeC4_ = 0;

    delete[] drudeC5_;
    drudeC5_ = 0;
  }
#else
  if (vcdt_)
  {
    delete[] vcdt_;
    vcdt_ = 0;

    delete[] omegasq_;
    omegasq_ = 0;

    delete[] eps_inf_;
    eps_inf_ = 0;
  }
#endif

  if (debyeA_)
  {
    delete[] debyeA_;
    debyeA_ = 0;
    delete[] debyeB_;
    debyeB_ = 0;
    delete[] debyeC_;
    debyeC_ = 0;
  }

  inited_ = false;
}

UPmlCommon::~UPmlCommon()
{
  deinit();
}
