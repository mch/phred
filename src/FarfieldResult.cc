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

#include "FarfieldResult.hh"
#include "Exceptions.hh"
#include "Constants.hh"
#include <math.h>

/* Jan's complex number and vector helper functions */
#include "fdtd_complex.hh"

FarfieldResult::FarfieldResult()
  : theta_start_(0), theta_stop_(0), phi_start_(90), phi_stop_(90),
    num_pts_(0), axis_(X_AXIS), freq_start_(0), freq_stop_(100), 
    num_freqs_(10), freq_space_(0), 
    theta_(0), phi_(0), 
    e_theta_re_(0), e_theta_im_(0), 
    e_phi_re_(0), e_phi_im_(0), 
    jff_mom_(0), mff_mom_(0),
    t_cross_(0), rank_(0), size_(0), 
    output_type_(ETHETAPHI), result_(0)
{
  freqs_.set_element_type(MPI_FLOAT);
  angles_.set_element_type(MPI_FLOAT);
  data_.set_element_type(MPI_FLOAT);

  freqs_.set_datatype(MPI_FLOAT);
  angles_.set_datatype(MPI_FLOAT);
  // data_'s datatype is different, calculated in init().

  data_.set_name("ff_data");

  angles_.has_time_dimension(false);
  freqs_.has_time_dimension(false);
  data_.has_time_dimension(false);

  variables_["angles"] = &angles_;
  variables_["freqs"] = &freqs_;
  variables_["ff_data"] = &data_;
}

FarfieldResult::~FarfieldResult()
{}

map<string, Variable *> & FarfieldResult::get_result(const Grid &grid, 
                                                     unsigned int time_step)
{
  if (result_time(time_step))
  {
    // Fetch the required data from other ranks. 

    // Do the computation
    ffpu_calculate(grid, time_step);

  }

  // Only return results on the last time step. 
  if (time_step == time_stop_)
  {
    // copy the data into the result_
    field_t *p = result_;
    for (unsigned int j = 0; j < num_freqs_; j++)
    {
      for (unsigned int i = 0; i < num_pts_; i++)
      {
        switch (output_type_)
        {
        case ETHETA:
          *p = e_theta_re_[i][j];
          ++p;
          *p = e_theta_im_[i][j];
          ++p;
          break;

        case EPHI:
          *p = e_phi_re_[i][j];
          ++p;
          *p = e_phi_im_[i][j];
          ++p;
          break;

        case HTHETA:
          *p = e_theta_re_[i][j] / ETA_0;
          ++p;
          *p = e_theta_im_[i][j] / ETA_0;
          ++p;
          break;

        case HPHI:
          *p = e_phi_re_[i][j] / ETA_0;
          ++p;
          *p = e_phi_im_[i][j] / ETA_0;
          ++p;
          break;

        case ETHETAPHI:
          *p = e_theta_re_[i][j];
          ++p;
          *p = e_theta_im_[i][j];
          ++p;
          *p = e_phi_re_[i][j];
          ++p;
          *p = e_phi_im_[i][j];
          ++p;
          break;

        case HTHETAPHI:
          *p = e_theta_re_[i][j] / ETA_0;
          ++p;
          *p = e_theta_im_[i][j] / ETA_0;
          ++p;
          *p = e_phi_re_[i][j] / ETA_0;
          ++p;
          *p = e_phi_im_[i][j] / ETA_0;
          ++p; 
          break;

        case RCSNORM:
          *p = freq_start_ + (j * freq_space_);
          ++p;
          *p = (e_theta_im_[i][j] * e_theta_im_[i][j] 
                + e_phi_im_[i][j] * e_phi_im_[i][j]
                + e_theta_re_[i][j] * e_theta_re_[i][j] 
                + e_phi_re_[i][j] * e_phi_re_[i][j]) 
            / (2 * ETA_0);
          ++p;
          break;

        case RCSDBPOL:
          *p = freq_start_ + (j * freq_space_);
          ++p;
          *p = 10 * log10((e_theta_im_[i][j] * e_theta_im_[i][j] 
                           + e_phi_im_[i][j] * e_phi_im_[i][j]
                           + e_theta_re_[i][j] * e_theta_re_[i][j] 
                           + e_phi_re_[i][j] * e_phi_re_[i][j]) 
                          / (2 * ETA_0));
          ++p;
          break;
        }
      }
    }

    data_.set_num(1);
    data_.set_ptr(result_);

    angles_.set_num(num_pts_);
    freqs_.set_num(num_freqs_);
  
    angles_.set_ptr(theta_);
    freqs_.set_ptr(freqs_buffer_);
  }

  return variables_;
}

void FarfieldResult::init(const Grid &grid)
{
  if (num_pts_ <= 0 || num_freqs_ <= 0)
    throw ResultException("A positive non-zero number of angle and frequency points are required.");

  if (size_ == 0)
    throw ResultException("The rank and size must be set in this result.");

  if (size_ > 1)
    throw ResultException("This result is not yet ready to be run in parallel.");

  t_cross_ = static_cast<int>(
    ceil(sqrt(pow(static_cast<delta_t>(grid.get_deltax() 
                                       * (global_r_.xmax - global_r_.xmin)),
                  static_cast<delta_t>(2.))
              + pow(static_cast<delta_t>(grid.get_deltay() 
                                         * (global_r_.ymax - global_r_.ymin)),
                    static_cast<delta_t>(2.))
              + pow(static_cast<delta_t>(grid.get_deltaz() 
                                         * (global_r_.zmax - global_r_.zmin)),
                    static_cast<delta_t>(2.)))
         / (C * grid.get_deltat())));
  
  theta_ = new float[num_pts_ * 2];
  phi_ = theta_ + num_pts_;
  freqs_buffer_ = new float[num_freqs_];

  e_theta_re_ = new float*[num_pts_];
  e_theta_im_ = new float*[num_pts_];
  e_phi_re_ = new float*[num_pts_];
  e_phi_im_ = new float*[num_pts_];
  
  jff_mom_ = new float**[num_pts_];
  mff_mom_ = new float**[num_pts_];
  
  if (!theta_ || !phi_ || !freqs_buffer_ || !e_theta_re_ 
      || !e_theta_im_ || !e_phi_re_
      || !e_phi_im_ || !jff_mom_ || !mff_mom_)
  {
    deinit(grid);
    throw MemoryException();
  }

  for (int i = 0; i < num_pts_; i++)
  {
    e_theta_re_[i] = new float[num_freqs_];
    e_theta_im_[i] = new float[num_freqs_];
    e_phi_re_[i] = new float[num_freqs_];
    e_phi_im_[i] = new float[num_freqs_];
    
    jff_mom_[i] = new float*[t_cross_];
    mff_mom_[i] = new float*[t_cross_];
    
    if (!e_theta_re_[i] || !e_theta_im_[i] || !e_phi_re_[i] || !e_phi_im_[i]
        || !jff_mom_[i] || !mff_mom_[i])
    {
      deinit(grid);
      throw MemoryException();
    }

    for (int j = 0; j < t_cross_; j++)
    {
      jff_mom_[i][j] = new float[3];
      mff_mom_[i][j] = new float[3];

      if (!jff_mom_[i][j] || !mff_mom_[i][j])
      {
        deinit(grid);
        throw MemoryException();
      }
    }
  }
  
  arc_connect();

  freq_space_ = (freq_stop_ - freq_start_) / num_freqs_;

  for (int i = 0; i < num_freqs_; i++)
    freqs_buffer_[i] = freq_start_ + i * freq_space_;

  // Output will be a three dimensional result, with a 2d table for
  // each frequency, and each frequency on it's own page. It will look
  // like this:
  // theta phi e_theta e_phi h_theta h_phi
  // or a variation on that. 
  unsigned int num_cols = 0;
  switch (output_type_)
  {
  case ETHETA:
    data_.set_name("E_theta");
    num_cols = 2;
    break;

  case EPHI:
    data_.set_name("E_phi");
    num_cols = 2;
    break;

  case HTHETA:
    data_.set_name("H_theta");
    num_cols = 2;
    break;

  case HPHI:
    data_.set_name("H_phi");
    num_cols = 2;
    break;

  case RCSNORM:
    data_.set_name("Norm_RCS");
    num_cols = 2;
    break;

  case RCSDBPOL:
    data_.set_name("Norm_RCS_db");
    num_cols = 2;
    break;

  case ETHETAPHI:
    data_.set_name("E_theta_phi");
    num_cols = 4;
    break;

  case HTHETAPHI:
    data_.set_name("H_theta_phi"); 
    num_cols = 4;
    break;
  }
  
  freqs_.add_dimension("frequency", num_freqs_, 0);

  angles_.add_dimension("theta and phi", 2, 0);
  angles_.add_dimension("angles", num_pts_, 0);

  data_.add_dimension("data", num_cols, 0);
  data_.add_dimension("points", num_pts_, 0);
  data_.add_dimension("frequencies", num_freqs_, 0);

  result_ = new field_t[num_cols * num_freqs_ * num_pts_];
  if (!result_)
  {
    deinit(grid);
    throw MemoryException();
  }
  memset(result_, 0, num_cols * num_freqs_ * num_pts_);
  
  MPI_Datatype temp;
  MPI_Type_contiguous(num_cols * num_freqs_ * num_pts_, 
                      MPI_FLOAT, &temp);
  MPI_Type_commit(&temp);
  data_.set_datatype(temp);
  data_.set_num(0);

  MPI_Datatype angle_type;
  MPI_Type_vector(num_pts_, 1, num_pts_, MPI_FLOAT, &angle_type);
  MPI_Type_commit(&angle_type);
  
}
  
void FarfieldResult::deinit(const Grid &grid)
{
  for (int i = 0; i < num_pts_; i++)
  {
    for (int j = 0; j < t_cross_; j++)
    {
      if (jff_mom_[i][j])
      {
        delete[] jff_mom_[i][j];
        jff_mom_[i][j] = 0;
      }

      if (mff_mom_[i][j])
      {
        delete[] mff_mom_[i][j];
        mff_mom_[i][j] = 0;
      }
    }

    if (e_theta_re_[i])
    {
      delete[] e_theta_re_[i];
      e_theta_re_[i] = 0;
    }

    if (e_theta_im_[i])
    {
      delete[] e_theta_im_[i];
      e_theta_im_[i] = 0;
    }

    if (e_phi_re_[i])
    {
      delete[] e_phi_re_[i];
      e_phi_re_[i] = 0;
    }

    if (e_phi_im_[i])
    {
      delete[] e_phi_im_[i];
      e_phi_im_[i] = 0;
    }

    if (jff_mom_[i])
    {
      delete[] jff_mom_[i];
      jff_mom_[i] = 0;
    }
    
    if (mff_mom_[i])
    {
      delete[] mff_mom_[i];    
      mff_mom_[i] = 0;    
    }
  }
  
  if (e_theta_re_)
  {
    delete[] e_theta_re_;
    e_theta_re_ = 0;
  }

  if (e_theta_im_)
  {
    delete[] e_theta_im_;
    e_theta_im_ = 0;
  }

  if (e_phi_re_)
  {
    delete[] e_phi_re_;
    e_phi_re_ = 0;
  }

  if (e_phi_im_)
  {
    delete[] e_phi_im_;
    e_phi_im_ = 0;
  }

  if (jff_mom_)
  {
    delete[] jff_mom_;
    jff_mom_ = 0;
  }
    
  if (mff_mom_)
  {
    delete[] mff_mom_;    
    mff_mom_ = 0;    
  }  

  if (theta_)
  {
    delete[] theta_;
    theta_ = 0;
  }

  if (phi_)
  {
    delete[] phi_;
    phi_ = 0;
  }

  if (result_)
  {
    delete[] result_;
    result_ = 0;
  }

  if (freqs_buffer_)
  {
    delete[] freqs_buffer_;
    freqs_buffer_ = 0;
  }
}

void FarfieldResult::arc_connect()
{
  // From Jan's implementation...
  float farf_delt;
  float farf_angle;
  float x_start,x_stop;
  float y_start,y_stop;
  float z_start,z_stop;
  int ptnr;

  if(num_pts_ == 1)
  {
    theta_[0] = theta_start_/180.*PI;
    phi_[0] = phi_start_/180.*PI;
  }
                
  else if(is_parallel( theta_start_, phi_start_,
                       theta_stop_, phi_stop_))
  {

    farf_delt = 2.*PI/( num_pts_-1);
    x_start = sin( theta_start_/180.*PI)*
      cos( phi_start_/180.*PI);      
    y_start = sin( theta_start_/180.*PI)*
      sin( phi_start_/180.*PI);      
    z_start = cos( theta_start_/180.*PI);
                                        
    for(ptnr = 0; ptnr < num_pts_; ptnr++)
    {
      switch (axis_)
      {
      case 1 :
        theta_[ptnr] = acos(
                           z_start*cos(ptnr*farf_delt)+
                           y_start*sin(ptnr*farf_delt));
        phi_[ptnr] = atan2_local(
                                y_start*cos(ptnr*farf_delt)-
                                z_start*sin(ptnr*farf_delt),x_start)
          ;
        break;
      case -1 :
        theta_[ptnr] = acos(
                           z_start*cos(-ptnr*farf_delt)+
                           y_start*sin(-ptnr*farf_delt));
        phi_[ptnr] = atan2_local(
                                y_start*cos(-ptnr*farf_delt)-
                                z_start*sin(-ptnr*farf_delt),x_start
                                );
        break;
      case 2 :
        theta_[ptnr] = acos(
                           x_start*cos(ptnr*farf_delt)+
                           z_start*sin(ptnr*farf_delt));
        phi_[ptnr] = atan2_local(
                                y_start,z_start*cos(ptnr*farf_delt)-
                                x_start*sin(ptnr*farf_delt));
        break;
      case -2 :
        theta_[ptnr] = acos(
                           x_start*cos(-ptnr*farf_delt)+
                           z_start*sin(-ptnr*farf_delt));
        phi_[ptnr] = atan2_local(
                                y_start,z_start*cos(-ptnr*farf_delt)
                                -
                                x_start*sin(-ptnr*farf_delt));
        break;
      case 3 :
        theta_[ptnr] = 
          theta_start_/180.*PI;
        phi_[ptnr] = 
          phi_start_/180.*PI
          +ptnr*farf_delt;
        break;
      default :
        theta_[ptnr] = 
          theta_start_/180.*PI;
        phi_[ptnr] = 
          phi_start_/180.*PI
          -ptnr*farf_delt;
        break;
      }
    }
  }
  else
  {
    x_start = sin( theta_start_/180.*PI)*
      cos( phi_start_/180.*PI);      
    y_start = sin( theta_start_/180.*PI)*
      sin( phi_start_/180.*PI);      
    z_start = cos( theta_start_/180.*PI);
                                                
    x_stop = sin( theta_stop_/180.*PI)*
      cos( phi_stop_/180.*PI);       
    y_stop = sin( theta_stop_/180.*PI)*
      sin( phi_stop_/180.*PI);       
    z_stop = cos( theta_stop_/180.*PI);
                                                
    farf_angle = acos(x_start*x_stop+y_start*y_stop+z_start*z_stop);
    if( axis_ < 0)
      farf_angle = 2.*PI-farf_angle;

    farf_delt = farf_angle/( num_pts_-1);
    for(ptnr=0;ptnr< num_pts_;ptnr++)
    {
      theta_[ptnr] = acos(1./sin(farf_angle)*
                         (sin(farf_angle-ptnr*farf_delt)*
                          z_start+sin(ptnr*farf_delt)*z_stop));
      phi_[ptnr] = atan2_local(
                              sin(farf_angle-ptnr*farf_delt)*y_start
                              +sin(ptnr*farf_delt)*y_stop,
                              sin(farf_angle-ptnr*farf_delt)*x_start
                              +sin(ptnr*farf_delt)*x_stop);
    }
  }
  return;
}

int FarfieldResult::is_parallel(float th1, float ph1, float th2, float ph2)
{
  int equal;
  int opposite;
  
  equal = (th1 == th2 && ph1 == ph2) || (th1 == 0 && th2 == 0) 
    || (th1 == 90 && th2 == 90);
  opposite = ((th1+th2) == 180 && fabs(ph1-ph2) == 180) 
    || (th1 == 180 && th2 == 0) || (th1 == 0 && th2 == 180);

  return(equal || opposite);
}

float FarfieldResult::atan2_local(float y, float x)
{
  if (x == 0 && y == 0)
    return(0.);
  else
    return(atan2(y, x));
}


void FarfieldResult::ffpu_calculate(const Grid &grid, 
                                    unsigned int time_step)
{
  int ptnr;
  float e_theta;
  float e_phi;
  float tt;
  int tnr,fnr;
  float omega;
  float tmp_im;
	
  ffpu_moment_update(grid);
	
  for(ptnr = 0; ptnr<num_pts_;ptnr++)
  {
    ffpu_ntff(theta_[ptnr],phi_[ptnr],
              mff_mom_[ptnr][0],jff_mom_[ptnr][0]
              ,&e_theta,&e_phi);
    
    if(num_freqs_ == 0)
    {
      e_theta_re_[ptnr][time_step-time_start_] =
        e_theta;
      e_phi_re_[ptnr][time_step-time_start_] =
        e_phi;
    }
    else
    {
      tt = time_step*grid.get_deltat();
      for(fnr = 0; fnr < num_freqs_ ; fnr++)
      {
        omega = 2*M_PI*(freq_start_
                        +fnr*freq_space_);
        e_theta_re_[ptnr][fnr] += 
          cos(-omega*tt)*e_theta;
        e_theta_im_[ptnr][fnr] += 
          sin(-omega*tt)*e_theta;
        e_phi_re_[ptnr][fnr] += 
          cos(-omega*tt)*e_phi;
        e_phi_im_[ptnr][fnr] += 
          sin(-omega*tt)*e_phi;
      }
    }
  }
		
  if (time_step == time_stop_) 
  {
    for(ptnr = 0; ptnr<num_pts_;ptnr++)
    {
      if(num_freqs_ == 0)
      {
        for(tnr = 0;tnr < time_stop_-
              time_start_+1; tnr++)
        
          e_theta_im_[ptnr][tnr] = (tnr == 0 ? 0. :
                                    (e_theta_re_[ptnr][tnr]
                                     -e_theta_re_[ptnr][tnr-1])
                                    /grid.get_deltat());
      
        for(tnr = 0;tnr < time_stop_-
              time_start_+1; tnr++)
        
          e_phi_im_[ptnr][tnr] = (tnr == 0 ? 0. :
                                  (e_phi_re_[ptnr][tnr]
                                   -e_phi_re_[ptnr][tnr-1])
                                  /grid.get_deltat());
      }
      else
        for(fnr = 0; fnr < num_freqs_ ; fnr++)
        {
          omega = 2*M_PI*(freq_start_
                          +fnr*freq_space_);
        
          tmp_im = 2.*sin(0.5*omega*grid.get_deltat())
            *e_theta_re_[ptnr][fnr];
        
          e_theta_re_[ptnr][fnr] = 
            -2.*sin(0.5*omega*grid.get_deltat())*
            e_theta_im_[ptnr][fnr];
        
          e_theta_im_[ptnr][fnr] = tmp_im;
	
          tmp_im = 2.*sin(0.5*omega*grid.get_deltat())
            *e_phi_re_[ptnr][fnr];
        
          e_phi_re_[ptnr][fnr] = 
            -2.*sin(0.5*omega*grid.get_deltat())*
            e_phi_im_[ptnr][fnr];
        
          e_phi_im_[ptnr][fnr] = tmp_im;
        }
    }
  }

  ffpu_moment_cycle();
	
}

void FarfieldResult::ffpu_moment_cycle()
{
  int ptnr,t;

	
  for(ptnr=0;ptnr<num_pts_;ptnr++)
    for(t=0;t<t_cross_-1;t++)
    {
      mff_mom_[ptnr][t][0] = mff_mom_[ptnr][t+1][0];
      mff_mom_[ptnr][t][1] = mff_mom_[ptnr][t+1][1];
      mff_mom_[ptnr][t][2] = mff_mom_[ptnr][t+1][2];
      jff_mom_[ptnr][t][0] = jff_mom_[ptnr][t+1][0];
      jff_mom_[ptnr][t][1] = jff_mom_[ptnr][t+1][1];
      jff_mom_[ptnr][t][2] = jff_mom_[ptnr][t+1][2];
    }
			
  for(ptnr=0;ptnr<num_pts_;ptnr++)
  {
    mff_mom_[ptnr][t_cross_-1][0] = 0.;
    mff_mom_[ptnr][t_cross_-1][1] = 0.;
    mff_mom_[ptnr][t_cross_-1][2] = 0.;
    jff_mom_[ptnr][t_cross_-1][0] = 0.;
    jff_mom_[ptnr][t_cross_-1][1] = 0.;
    jff_mom_[ptnr][t_cross_-1][2] = 0.;
  }
  return;
}

void FarfieldResult::ffpu_ntff(float theta, float phi, float *ms, float *js,
                               float *e_theta, float *e_phi)
{
  float theta_unit[3];
  float phi_unit[3];
  float ms_theta,ms_phi;
  float js_theta,js_phi;
	
	
  f_vfill(theta_unit,cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta));
  f_vfill(phi_unit,-sin(phi),cos(phi),0.);
	
  ms_theta = f_dot(ms,theta_unit);
  ms_phi = f_dot(ms,phi_unit);
	
  js_theta = f_dot(js,theta_unit);
  js_phi = f_dot(js,phi_unit);
	
  *e_theta = (-ms_phi-ETA_0*js_theta)/(4.*PI*C);
  *e_phi = (ms_theta-ETA_0*js_phi)/(4.*PI*C);
  return;
}

void FarfieldResult::ffpu_moment_update(const Grid &grid)
{
  int surf_i,surf_j,surf_k;
  int glob_i,glob_j,glob_k;
  float r_unit[3];
  float r_accent[3];
  float origin[3];
  float e_x,e_y,e_z;
  float h_x,h_y,h_z;
  int ptnr;
  float tau;
  int t1_e,t1_h;
  float f1_e,f1_h;
  float f2_e,f2_h;
	
  for(ptnr = 0; ptnr < num_pts_; ptnr++)
  {
    r_unit[0] = sin(theta_[ptnr])*cos(phi_[ptnr]);
    r_unit[1] = sin(theta_[ptnr])*sin(phi_[ptnr]);
    r_unit[2] = cos(theta_[ptnr]);
    origin[0] = (r_unit[0] < 0. ? (global_r_.xmin-1)*grid.get_deltax() :
                 (global_r_.xmax+1)*grid.get_deltax());
    origin[1] = (r_unit[1] < 0. ? (global_r_.ymin-1)*grid.get_deltay() :
                 (global_r_.ymax+1)*grid.get_deltay());
    origin[2] = (r_unit[2] < 0. ? (global_r_.zmin-1)*grid.get_deltaz() :
                 (global_r_.zmax+1)*grid.get_deltaz());
	

    /* xmin side */
	
    for(surf_j=0, glob_j=global_r_.ymin+surf_j ;
        surf_j<global_r_.ymax-global_r_.ymin ;
        surf_j++,glob_j++)
      for(surf_k=0, glob_k=global_r_.zmin+surf_k ;
          surf_k<global_r_.zmax-global_r_.zmin ;
          surf_k++,glob_k++)
      {
        r_accent[0] = global_r_.xmin*grid.get_deltax()-origin[0];
        r_accent[1] = (glob_j+0.5)*grid.get_deltay()-origin[1];
        r_accent[2] = (glob_k+0.5)*grid.get_deltaz()-origin[2];
        tau = -f_dot(r_unit,r_accent)/(C*grid.get_deltat());

        t1_h = static_cast<int>(floor(tau));
        f2_h = tau-t1_h;
        f1_h = 1.-f2_h;

        t1_e = static_cast<int>(floor(tau+0.5));
        f2_e = tau+0.5-t1_e;
        f1_e = 1.-f2_e;
				
        
        e_y = 0.5*(grid.get_ey(global_r_.xmin,glob_j,glob_k)+
                   grid.get_ey(global_r_.xmin,glob_j,glob_k+1));
        e_z = 0.5*(grid.get_ez(global_r_.xmin,glob_j,glob_k)+
                   grid.get_ez(global_r_.xmin,glob_j+1,glob_k));
        h_y = 0.25*(grid.get_hy(global_r_.xmin-1,glob_j,glob_k)+
                    grid.get_hy(global_r_.xmin-1,glob_j+1,glob_k)+
                    grid.get_hy(global_r_.xmin,glob_j,glob_k)+
                    grid.get_hy(global_r_.xmin,glob_j+1,glob_k));
        h_z = 0.25*(grid.get_hz(global_r_.xmin-1,glob_j,glob_k)+
                    grid.get_hz(global_r_.xmin-1,glob_j,glob_k+1)+
                    grid.get_hz(global_r_.xmin,glob_j,glob_k)+
                    grid.get_hz(global_r_.xmin,glob_j,glob_k+1));

        mff_mom_[ptnr][t1_e][1] -= f1_e*e_z*grid.get_deltay()*grid.get_deltaz();
        mff_mom_[ptnr][t1_e+1][1] -= f2_e*e_z*grid.get_deltay()*grid.get_deltaz();
        mff_mom_[ptnr][t1_e][2] += f1_e*e_y*grid.get_deltay()*grid.get_deltaz();
        mff_mom_[ptnr][t1_e+1][2] += f2_e*e_y*grid.get_deltay()*grid.get_deltaz();
        jff_mom_[ptnr][t1_h][1] += f1_h*h_z*grid.get_deltay()*grid.get_deltaz();
        jff_mom_[ptnr][t1_h+1][1] += f2_h*h_z*grid.get_deltay()*grid.get_deltaz();
        jff_mom_[ptnr][t1_h][2] -= f1_h*h_y*grid.get_deltay()*grid.get_deltaz();
        jff_mom_[ptnr][t1_h+1][2] -= f2_h*h_y*grid.get_deltay()*grid.get_deltaz();
      }
	
    /* xmax side */
	
    for(surf_j=0, glob_j = global_r_.ymin + surf_j;
        surf_j < global_r_.ymax - global_r_.ymin;
        surf_j++, glob_j++)
      for(surf_k = 0, glob_k = global_r_.zmin + surf_k;
          surf_k < global_r_.zmax - global_r_.zmin;
          surf_k++,glob_k++)
      {
        r_accent[0] = global_r_.xmax*grid.get_deltax()-origin[0];
        r_accent[1] = (glob_j+0.5)*grid.get_deltay()-origin[1];
        r_accent[2] = (glob_k+0.5)*grid.get_deltaz()-origin[2];
        tau = -(r_unit[0]*r_accent[0]+r_unit[1]*r_accent[1]+
                r_unit[2]*r_accent[2])/(C*grid.get_deltat());

        t1_h = static_cast<int>(floor(tau));
        f2_h = tau-t1_h;
        f1_h = 1.-f2_h;

        t1_e = static_cast<int>(floor(tau+0.5));
        f2_e = tau+0.5-t1_e;
        f1_e = 1.-f2_e;
				
        e_y = 0.5*(grid.get_ey(global_r_.xmax,glob_j,glob_k)+
                   grid.get_ey(global_r_.xmax,glob_j,glob_k+1));
        e_z = 0.5*(grid.get_ez(global_r_.xmax,glob_j,glob_k)+
                   grid.get_ez(global_r_.xmax,glob_j+1,glob_k));
        h_y = 0.25*(grid.get_hy(global_r_.xmax-1,glob_j,glob_k)+
                    grid.get_hy(global_r_.xmax-1,glob_j+1,glob_k)+
                    grid.get_hy(global_r_.xmax,glob_j,glob_k)+
                    grid.get_hy(global_r_.xmax,glob_j+1,glob_k));
        h_z = 0.25*(grid.get_hz(global_r_.xmax-1,glob_j,glob_k)+
                    grid.get_hz(global_r_.xmax-1,glob_j,glob_k+1)+
                    grid.get_hz(global_r_.xmax,glob_j,glob_k)+
                    grid.get_hz(global_r_.xmax,glob_j,glob_k+1));

        mff_mom_[ptnr][t1_e][1] += f1_e*e_z*grid.get_deltay()*grid.get_deltaz();
        mff_mom_[ptnr][t1_e+1][1] += f2_e*e_z*grid.get_deltay()*grid.get_deltaz();
        mff_mom_[ptnr][t1_e][2] -= f1_e*e_y*grid.get_deltay()*grid.get_deltaz();
        mff_mom_[ptnr][t1_e+1][2] -= f2_e*e_y*grid.get_deltay()*grid.get_deltaz();
        jff_mom_[ptnr][t1_h][1] -= f1_h*h_z*grid.get_deltay()*grid.get_deltaz();
        jff_mom_[ptnr][t1_h+1][1] -= f2_h*h_z*grid.get_deltay()*grid.get_deltaz();
        jff_mom_[ptnr][t1_h][2] += f1_h*h_y*grid.get_deltay()*grid.get_deltaz();
        jff_mom_[ptnr][t1_h+1][2] += f2_h*h_y*grid.get_deltay()*grid.get_deltaz();
      }
		
    /* side ymin */	
	
    for(surf_i=0, glob_i=global_r_.xmin+surf_i ;
        surf_i<global_r_.xmax-global_r_.xmin ;
        surf_i++,glob_i++)
      for(surf_k=0, glob_k=global_r_.zmin+surf_k ;
          surf_k<global_r_.zmax-global_r_.zmin ;
          surf_k++,glob_k++)
      {
        r_accent[0] = (glob_i+0.5)*grid.get_deltax()-origin[0];
        r_accent[1] = global_r_.ymin*grid.get_deltay()-origin[1];
        r_accent[2] = (glob_k+0.5)*grid.get_deltaz()-origin[2];

        tau = -(r_unit[0]*r_accent[0]+r_unit[1]*r_accent[1]+
                r_unit[2]*r_accent[2])/(C*grid.get_deltat());

        t1_h = static_cast<int>(floor(tau));
        f2_h = tau-t1_h;
        f1_h = 1.-f2_h;

        t1_e = static_cast<int>(floor(tau+0.5));
        f2_e = tau+0.5-t1_e;
        f1_e = 1.-f2_e;
				
        e_x = 0.5*(grid.get_ex(glob_i,global_r_.ymin,glob_k)+
                   grid.get_ex(glob_i,global_r_.ymin,glob_k+1));
        e_z = 0.5*(grid.get_ez(glob_i,global_r_.ymin,glob_k)+
                   grid.get_ez(glob_i+1,global_r_.ymin,glob_k));
        h_x = 0.25*(grid.get_hx(glob_i,global_r_.ymin-1,glob_k)+
                    grid.get_hx(glob_i+1,global_r_.ymin-1,glob_k)+
                    grid.get_hx(glob_i,global_r_.ymin,glob_k)+
                    grid.get_hx(glob_i+1,global_r_.ymin,glob_k));
        h_z = 0.25*(grid.get_hz(glob_i,global_r_.ymin-1,glob_k)+
                    grid.get_hz(glob_i,global_r_.ymin-1,glob_k+1)+
                    grid.get_hz(glob_i,global_r_.ymin,glob_k)+
                    grid.get_hz(glob_i,global_r_.ymin,glob_k+1));
	
        mff_mom_[ptnr][t1_e][0] += f1_e*e_z*grid.get_deltax()*grid.get_deltaz();
        mff_mom_[ptnr][t1_e+1][0] += f2_e*e_z*grid.get_deltax()*grid.get_deltaz();
        mff_mom_[ptnr][t1_e][2] -= f1_e*e_x*grid.get_deltax()*grid.get_deltaz();
        mff_mom_[ptnr][t1_e+1][2] -= f2_e*e_x*grid.get_deltax()*grid.get_deltaz();
        jff_mom_[ptnr][t1_h][0] -= f1_h*h_z*grid.get_deltax()*grid.get_deltaz();
        jff_mom_[ptnr][t1_h+1][0] -= f2_h*h_z*grid.get_deltax()*grid.get_deltaz();
        jff_mom_[ptnr][t1_h][2] += f1_h*h_x*grid.get_deltax()*grid.get_deltaz();
        jff_mom_[ptnr][t1_h+1][2] += f2_h*h_x*grid.get_deltax()*grid.get_deltaz();
      }
				
    /* side ymax */	
	
    for(surf_i=0, glob_i=global_r_.xmin+surf_i ;
        surf_i<global_r_.xmax-global_r_.xmin ;
        surf_i++,glob_i++)
      for(surf_k=0, glob_k=global_r_.zmin+surf_k ;
          surf_k<global_r_.zmax-global_r_.zmin ;
          surf_k++,glob_k++)
      {
        r_accent[0] = (glob_i+0.5)*grid.get_deltax()-origin[0];
        r_accent[1] = global_r_.ymax*grid.get_deltay()-origin[1];
        r_accent[2] = (glob_k+0.5)*grid.get_deltaz()-origin[2];

        tau = -(r_unit[0]*r_accent[0]+r_unit[1]*r_accent[1]+
                r_unit[2]*r_accent[2])/(C*grid.get_deltat());

        t1_h = static_cast<int>(floor(tau));
        f2_h = tau-t1_h;
        f1_h = 1.-f2_h;

        t1_e = static_cast<int>(floor(tau+0.5));
        f2_e = tau+0.5-t1_e;
        f1_e = 1.-f2_e;
				
        e_x = 0.5*(grid.get_ex(glob_i,global_r_.ymax,glob_k)+
                   grid.get_ex(glob_i,global_r_.ymax,glob_k+1));
        e_z = 0.5*(grid.get_ez(glob_i,global_r_.ymax,glob_k)+
                   grid.get_ez(glob_i+1,global_r_.ymax,glob_k));
        h_x = 0.25*(grid.get_hx(glob_i,global_r_.ymax-1,glob_k)+
                    grid.get_hx(glob_i+1,global_r_.ymax-1,glob_k)+
                    grid.get_hx(glob_i,global_r_.ymax,glob_k)+
                    grid.get_hx(glob_i+1,global_r_.ymax,glob_k));
        h_z = 0.25*(grid.get_hz(glob_i,global_r_.ymax-1,glob_k)+
                    grid.get_hz(glob_i,global_r_.ymax-1,glob_k+1)+
                    grid.get_hz(glob_i,global_r_.ymax,glob_k)+
                    grid.get_hz(glob_i,global_r_.ymax,glob_k+1));
	
        mff_mom_[ptnr][t1_e][0] -= f1_e*e_z*grid.get_deltax()*grid.get_deltaz();
        mff_mom_[ptnr][t1_e+1][0] -= f2_e*e_z*grid.get_deltax()*grid.get_deltaz();
        mff_mom_[ptnr][t1_e][2] += f1_e*e_x*grid.get_deltax()*grid.get_deltaz();
        mff_mom_[ptnr][t1_e+1][2] += f2_e*e_x*grid.get_deltax()*grid.get_deltaz();
        jff_mom_[ptnr][t1_h][0] += f1_h*h_z*grid.get_deltax()*grid.get_deltaz();
        jff_mom_[ptnr][t1_h+1][0] += f2_h*h_z*grid.get_deltax()*grid.get_deltaz();
        jff_mom_[ptnr][t1_h][2] -= f1_h*h_x*grid.get_deltax()*grid.get_deltaz();
        jff_mom_[ptnr][t1_h+1][2] -= f2_h*h_x*grid.get_deltax()*grid.get_deltaz();
      }
			
    /* side zmin */	
	
    for(surf_i=0, glob_i=global_r_.xmin+surf_i ;
        surf_i<global_r_.xmax-global_r_.xmin ;
        surf_i++,glob_i++)
      for(surf_j=0, glob_j=global_r_.ymin+surf_j ;
          surf_j<global_r_.ymax-global_r_.ymin ;
          surf_j++,glob_j++)
      {
        r_accent[0] = (glob_i+0.5)*grid.get_deltax()-origin[0];
        r_accent[1] = (glob_j+0.5)*grid.get_deltay()-origin[1];
        r_accent[2] = global_r_.zmin*grid.get_deltaz()-origin[2];

        tau = -(r_unit[0]*r_accent[0]+r_unit[1]*r_accent[1]+
                r_unit[2]*r_accent[2])/(C*grid.get_deltat());

        t1_h = static_cast<int>(floor(tau));
        f2_h = tau-t1_h;
        f1_h = 1.-f2_h;

        t1_e = static_cast<int>(floor(tau+0.5));
        f2_e = tau+0.5-t1_e;
        f1_e = 1.-f2_e;
				

        e_x = 0.5*(grid.get_ex(glob_i,glob_j,global_r_.zmin)+
                   grid.get_ex(glob_i,glob_j+1,global_r_.zmin));
        e_y = 0.5*(grid.get_ey(glob_i,glob_j,global_r_.zmin)+
                   grid.get_ey(glob_i+1,glob_j,global_r_.zmin));
        h_x = 0.25*(grid.get_hx(glob_i,glob_j,global_r_.zmin-1)+
                    grid.get_hx(glob_i+1,glob_j,global_r_.zmin-1)+
                    grid.get_hx(glob_i,glob_j,global_r_.zmin)+
                    grid.get_hx(glob_i+1,glob_j,global_r_.zmin));
        h_y = 0.25*(grid.get_hy(glob_i,glob_j,global_r_.zmin-1)+
                    grid.get_hy(glob_i,glob_j+1,global_r_.zmin-1)+
                    grid.get_hy(glob_i,glob_j,global_r_.zmin)+
                    grid.get_hy(glob_i,glob_j+1,global_r_.zmin));
        mff_mom_[ptnr][t1_e][0] -= f1_e*e_y*grid.get_deltax()*grid.get_deltay();
        mff_mom_[ptnr][t1_e+1][0] -= f2_e*e_y*grid.get_deltax()*grid.get_deltay();
        mff_mom_[ptnr][t1_e][1] += f1_e*e_x*grid.get_deltax()*grid.get_deltay();
        mff_mom_[ptnr][t1_e+1][1] += f2_e*e_x*grid.get_deltax()*grid.get_deltay();
        jff_mom_[ptnr][t1_h][0] += f1_h*h_y*grid.get_deltax()*grid.get_deltay();
        jff_mom_[ptnr][t1_h+1][0] += f2_h*h_y*grid.get_deltax()*grid.get_deltay();
        jff_mom_[ptnr][t1_h][1] -= f1_h*h_x*grid.get_deltax()*grid.get_deltay();
        jff_mom_[ptnr][t1_h+1][1] -= f2_h*h_x*grid.get_deltax()*grid.get_deltay();
      }			
		
    /* side zmax */	
	
    for(surf_i=0, glob_i=global_r_.xmin+surf_i ;
        surf_i<global_r_.xmax-global_r_.xmin ;
        surf_i++,glob_i++)
      for(surf_j=0, glob_j=global_r_.ymin+surf_j ;
          surf_j<global_r_.ymax-global_r_.ymin ;
          surf_j++,glob_j++)
      {
        r_accent[0] = (glob_i+0.5)*grid.get_deltax()-origin[0];
        r_accent[1] = (glob_j+0.5)*grid.get_deltay()-origin[1];
        r_accent[2] = global_r_.zmax*grid.get_deltaz()-origin[2];

        tau = -(r_unit[0]*r_accent[0]+r_unit[1]*r_accent[1]+
                r_unit[2]*r_accent[2])/(C*grid.get_deltat());

        t1_h = static_cast<int>(floor(tau));
        f2_h = tau-t1_h;
        f1_h = 1.-f2_h;

        t1_e = static_cast<int>(floor(tau+0.5));
        f2_e = tau+0.5-t1_e;
        f1_e = 1.-f2_e;
				

        e_x = 0.5*(grid.get_ex(glob_i,glob_j,global_r_.zmax)+
                   grid.get_ex(glob_i,glob_j+1,global_r_.zmax));
        e_y = 0.5*(grid.get_ey(glob_i,glob_j,global_r_.zmax)+
                   grid.get_ey(glob_i+1,glob_j,global_r_.zmax));
        h_x = 0.25*(grid.get_hx(glob_i,glob_j,global_r_.zmax-1)+
                    grid.get_hx(glob_i+1,glob_j,global_r_.zmax-1)+
                    grid.get_hx(glob_i,glob_j,global_r_.zmax)+
                    grid.get_hx(glob_i+1,glob_j,global_r_.zmax));
        h_y = 0.25*(grid.get_hy(glob_i,glob_j,global_r_.zmax-1)+
                    grid.get_hy(glob_i,glob_j+1,global_r_.zmax-1)+
                    grid.get_hy(glob_i,glob_j,global_r_.zmax)+
                    grid.get_hy(glob_i,glob_j+1,global_r_.zmax));
        mff_mom_[ptnr][t1_e][0] += f1_e*e_y*grid.get_deltax()*grid.get_deltay();
        mff_mom_[ptnr][t1_e+1][0] += f2_e*e_y*grid.get_deltax()*grid.get_deltay();
        mff_mom_[ptnr][t1_e][1] -= f1_e*e_x*grid.get_deltax()*grid.get_deltay();
        mff_mom_[ptnr][t1_e+1][1] -= f2_e*e_x*grid.get_deltax()*grid.get_deltay();
        jff_mom_[ptnr][t1_h][0] -= f1_h*h_y*grid.get_deltax()*grid.get_deltay();
        jff_mom_[ptnr][t1_h+1][0] -= f2_h*h_y*grid.get_deltax()*grid.get_deltay();
        jff_mom_[ptnr][t1_h][1] += f1_h*h_x*grid.get_deltax()*grid.get_deltay();
        jff_mom_[ptnr][t1_h+1][1] += f2_h*h_x*grid.get_deltax()*grid.get_deltay();
      }
  }
  return;
}
