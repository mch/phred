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
#include <math.h>

FarfieldResult::FarfieldResult()
  : theta_start_(0), theta_stop_(0), phi_start_(90), phi_stop_(90),
    num_pts_(0), axis_(X_AXIS), freq_start_(0), freq_stop_(100), 
    num_freqs_(10), freq_space_(0), 
    theta_(0), phi_(0), 
    e_theta_re_(0), e_theta_im_(0), 
    e_phi_re_(0), e_phi_im_(0), 
    jff_mom_(0), mff_mom_(0),
    t_cross_(0), rank_(0), size_(0)
{}

FarfieldResult~FarfieldResult()
{}

Data & FarfieldResult::get_result(const Grid &grid, unsigned int time_step)
{
  if (result_time(time_step))
  {
    // Fetch the required data from other ranks. 

    // Do the computation

    // Return the results...
    
  }
}

void FarfieldResult::init(const Grid &grid)
{
  if (num_pts_ <= 0 || num_freqs_ <= 0)
    throw ResultException("A positive non-zero number of angle and frequency points are required.");

  if (size_ == 0)
    throw ResultException("The rank and size must be set in this result.");

  if (size_ > 1)
    throw ResultException("This result is not yet ready to be run in parallel.");

  t_cross_ = 
    ceil(sqrt(pow(grid.get_deltax() * (global_r_.xmax - global_r_.xmin),2.)
              + pow(grid.get_deltay() * (global_r_.ymax - global_r_.ymin),2.)
              + pow(gridptr->deltaz * (global_r_.zmax - global_r_.zmin),2.))
         / (P_C * grid.get_deltat()));

  theta_ = new float[num_pts_];
  phi_ = new float[num_pts_];
  
  e_theta_re_ = new float*[num_pts_];
  e_theta_im_ = new float*[num_pts_];
  e_phi_re_ = new float*[num_pts_];
  e_phi_im_ = new float*[num_pts_];
  
  jff_mom_ = new float**[num_pts];
  mff_mom_ = new float**[num_pts];
  
  if (!theta_ || !phi_ || !e_theta_re_ || !e_theta_im_ || !e_phi_re_
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
    theta_[0] = theta_start/180.*PI;
    phi_[0] = phi_start/180.*PI;
  }
                
  else if(is_parallel( theta_start, phi_start,
                       theta_stop, phi_stop))
  {

    farf_delt = 2.*PI/( num_pts_-1);
    x_start = sin( theta_start/180.*PI)*
      cos( phi_start/180.*PI);      
    y_start = sin( theta_start/180.*PI)*
      sin( phi_start/180.*PI);      
    z_start = cos( theta_start/180.*PI);
                                        
    for(ptnr=0;ptnr< num_pts_;ptnr++)
    {
      switch( axis_of_rotation)
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
          theta_start/180.*PI;
        phi_[ptnr] = 
          phi_start/180.*PI
          +ptnr*farf_delt;
        break;
      default :
        theta_[ptnr] = 
          theta_start/180.*PI;
        phi_[ptnr] = 
          phi_start/180.*PI
          -ptnr*farf_delt;
        break;
      }
    }
  }
  else
  {
    x_start = sin( theta_start/180.*PI)*
      cos( phi_start/180.*PI);      
    y_start = sin( theta_start/180.*PI)*
      sin( phi_start/180.*PI);      
    z_start = cos( theta_start/180.*PI);
                                                
    x_stop = sin( theta_stop/180.*PI)*
      cos( phi_stop/180.*PI);       
    y_stop = sin( theta_stop/180.*PI)*
      sin( phi_stop/180.*PI);       
    z_stop = cos( theta_stop/180.*PI);
                                                
    farf_angle = acos(x_start*x_stop+y_start*y_stop+z_start*z_stop);
    if( axis_of_rotation < 0)
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
