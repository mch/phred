/* 
   Phred - Phred is a parallel finite difference time domain
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
#include "../config.h"
#include "../Exceptions.hh"
#include "../Constants.hh"
#include "../Globals.hh"
#include "../PlaneTiling.hh"
#include <cmath>
#include <cstring>
#include <cassert>

#include <fstream>

using namespace std;

void FarfieldResult::idx_tests()
{
  ofstream wuf;

  wuf.open("wutest.txt", ofstream::out);
  wuf << "# phis: " << phi_data_.length() << "\n";
  wuf << "# thetas: " << theta_data_.length() << "\n";
  wuf << "# fftsteps: " << ff_tsteps_ << "\n";

  wuf << "phi_idx phi theta_idx theta WU_index\n";

  for (int phi_idx = 0; phi_idx < phi_data_.length(); phi_idx++)
  {
    for (int theta_idx = 0; theta_idx < theta_data_.length(); theta_idx++)
    {
      wuf << phi_idx << "\t" << phi_data_.get(phi_idx) << "\t"
          << theta_idx << "\t" << theta_data_.get(theta_idx) << "\t"
          << WU_index(phi_idx, theta_idx, 0, theta_data_.length(),
                      ff_tsteps_) << "\n";
    }
  }

  wuf.close();
}

void FarfieldResult::dump_temps()
{
  ofstream fp;

  fp.open("cartesian_fields.txt", ofstream::out);
  
  int idx = 0;
  int t_start = (ff_tsteps_ / 2) - (ff_tsteps_ / 8);
  int t_end = (ff_tsteps_ / 2) + (ff_tsteps_ / 8);

  fp << scientific;

  for (int phi_idx = 0; phi_idx < phi_data_.length(); phi_idx++)
  {
    for (int theta_idx = 0; theta_idx < theta_data_.length(); theta_idx++)
    {
      for (int t = t_start; t < t_end; t++)
      {
        fp << phi_data_.get(phi_idx) << " " << theta_data_.get(theta_idx)
           << " " << t << " " 
           << Wx_[idx] << " " << Wy_[idx] << " " << Wz_[idx] << " "
           << Ux_[idx] << " " << Uy_[idx] << " " << Uz_[idx]
           << "\n";

        idx++;
      }
    }
  }

  fp.close();

}

FarfieldResult::FarfieldResult()
  : r_(100),
    Wx_(0), Wy_(0), Wz_(0), Ux_(0), Uy_(0), Uz_(0),
    E_theta_(0), E_phi_(0)
{
  for (int i = 0; i < 6; i++)
    use_face_[i] = true;
}

FarfieldResult::~FarfieldResult()
{
  deinit();
}

void FarfieldResult::set_radius(field_t r)
{
  if (r > 0.0)
  {
    r_ = r;
  } 
  else 
  {
    throw ResultException("Farfield radius must be greater than zero.");
  }
}

void FarfieldResult::set_theta_degrees(field_t theta_start, 
                                        field_t theta_stop, 
                                        unsigned int num_theta)
{
  set_theta(theta_start * (PI / 180), theta_stop * (PI / 180),
            num_theta);
}

void FarfieldResult::set_theta(field_t theta_start, field_t theta_stop, 
                                unsigned int num_theta)
{
  if (num_theta < 1)
    throw ResultException("FarfieldResult must return data at one "
                          "or more angles of theta");

  if (abs(theta_stop - theta_start) > 2 * PI)
    throw ResultException("FarfiedlResult theta angles must not span more "
                          "than 360 degrees.");

  theta_data_.set_params(theta_start, theta_stop, num_theta);
}

void FarfieldResult::set_phi_degrees(field_t phi_start, field_t phi_stop, 
                                      unsigned int num_phi)
{
  set_phi(phi_start * (PI / 180), phi_stop * (PI / 180),
          num_phi);
}

void FarfieldResult::set_phi(field_t phi_start, field_t phi_stop, 
                                unsigned int num_phi)
{
  if (num_phi < 1)
    throw ResultException("FarfieldResult must return data at one "
                          "or more angles of phi");
  
  if (abs(phi_stop - phi_start) > 2 * PI)
    throw ResultException("FarfiedlResult phi angles must not span more "
                          "than 360 degrees.");

  phi_data_.set_params(phi_start, phi_stop, num_phi);
}

void FarfieldResult::init(const Grid &grid)
{
  deinit();

  if (box_.get())
  {
    region_ = grid.get_local_region(*(box_.get()));
  } else {
    throw ResultException("SurfaceCurrentResult has no surface defined!");
  }

  shared_ptr<Block> gregion = grid.get_global_region(*(box_.get()));

  // Number of farfield timesteps we have to keep track of
  point bsz = (*box_).get_size();

  // Maximum distance between any two points on the Huygen's surface.
  field_t maxlength = sqrt(bsz.x * bsz.x + bsz.y * bsz.y + bsz.z * bsz.z);

  t_cross_ = maxlength / C;

  ff_tsteps_ = (time_stop_ - time_start_) 
    + 2 * static_cast<unsigned int>(ceil(t_cross_ / grid.get_deltat()));
  
#ifdef DEBUG
  cerr << "Farfield, num ff_tsteps: " << ff_tsteps_ << endl;
#endif

  // Allocate space for W and U
  unsigned int sz = ff_tsteps_ * theta_data_.length() * phi_data_.length();

  unsigned int dftsz = frequencies_.length() 
    * theta_data_.length() * phi_data_.length();

  Wx_ = new field_t[sz];
  Wy_ = new field_t[sz];
  Wz_ = new field_t[sz];
  
  Ux_ = new field_t[sz];
  Uy_ = new field_t[sz];
  Uz_ = new field_t[sz];

  memset(Wx_, 0, sizeof(field_t)*sz);
  memset(Wy_, 0, sizeof(field_t)*sz);
  memset(Wz_, 0, sizeof(field_t)*sz);

  memset(Ux_, 0, sizeof(field_t)*sz);
  memset(Uy_, 0, sizeof(field_t)*sz);
  memset(Uz_, 0, sizeof(field_t)*sz);

  if (MPI_RANK == 0)
  {
    E_theta_ = new field_t[sz];
    E_phi_ = new field_t[sz];

    memset(E_theta_, 0, sizeof(field_t)*sz);
    memset(E_phi_, 0, sizeof(field_t)*sz);

    e_theta_real_ = new field_t[dftsz];
    e_theta_imag_ = new field_t[dftsz];
    e_phi_real_ = new field_t[dftsz];
    e_phi_imag_ = new field_t[dftsz];
    rcs_ = new field_t[dftsz];

    memset(e_theta_real_, 0, sizeof(field_t)*dftsz);
    memset(e_theta_imag_, 0, sizeof(field_t)*dftsz);
    memset(e_phi_real_, 0, sizeof(field_t)*dftsz);
    memset(e_phi_imag_, 0, sizeof(field_t)*dftsz);
    memset(rcs_, 0, sizeof(field_t)*dftsz);
  }

  // Set up output variables
  freqs_.reset();
  theta_.reset();
  phi_.reset();
  E_theta_var_.reset();
  E_phi_var_.reset();
  e_tr_.reset();
  e_ti_.reset();
  e_pr_.reset();
  e_pi_.reset();
  rcs_var_.reset();

  pre_vars_["freqs"] = &freqs_;
  pre_vars_["theta"] = &theta_;
  pre_vars_["phi"] = &phi_;

  post_vars_["Et"] = &E_theta_var_;
  post_vars_["Ep"] = &E_phi_var_;
  post_vars_["Et_real"] = &e_tr_;
  post_vars_["Et_imag"] = &e_ti_;
  post_vars_["Ep_real"] = &e_pr_;
  post_vars_["Ep_imag"] = &e_pi_;
  post_vars_["rcs"] = &rcs_var_;

  // Frequency
  freqs_.add_dimension("freqs", frequencies_.length(), 
                       frequencies_.length(), 0);
  freqs_.set_name(base_name_ + "freqs");
  freqs_.set_ptr(frequencies_.get_ptr());
  freqs_.has_time_dimension(false);

  // Theta
  theta_.add_dimension("theta", theta_data_.length(), 
                       theta_data_.length(), 0);
  theta_.set_name(base_name_ + "theta");
  theta_.set_ptr(theta_data_.get_ptr());
  theta_.has_time_dimension(false);

  // Phi
  phi_.add_dimension("phi", phi_data_.length(), phi_data_.length(), 0);
  phi_.set_name(base_name_ + "phi");
  phi_.set_ptr(phi_data_.get_ptr());
  phi_.has_time_dimension(false);

  // Normalized Electric field in spherical coords
  E_theta_var_.add_dimension("phi", phi_data_.length(), 
                             phi_data_.length(), 0);
  E_theta_var_.add_dimension("theta", theta_data_.length(), 
                             theta_data_.length(), 0);
  E_theta_var_.add_dimension("fftime", ff_tsteps_, ff_tsteps_, 0);
  E_theta_var_.set_name(base_name_ + "Et");
  E_theta_var_.set_ptr(E_theta_);
  E_theta_var_.has_time_dimension(false);

  E_phi_var_.add_dimension("phi", phi_data_.length(), 
                             phi_data_.length(), 0);
  E_phi_var_.add_dimension("theta", theta_data_.length(), 
                             theta_data_.length(), 0);
  E_phi_var_.add_dimension("fftime", ff_tsteps_, ff_tsteps_, 0);
  E_phi_var_.set_name(base_name_ + "Ep");
  E_phi_var_.set_ptr(E_phi_);
  E_phi_var_.has_time_dimension(false);

  // DFT'd electric field
  e_tr_.add_dimension("phi", phi_data_.length(), 
                      phi_data_.length(), 0);
  e_tr_.add_dimension("theta", theta_data_.length(), 
                      theta_data_.length(), 0);
  e_tr_.add_dimension("freq", frequencies_.length(),
                      frequencies_.length(), 0);
  e_tr_.set_name(base_name_ + "Et_real");
  e_tr_.set_ptr(e_theta_real_);
  e_tr_.has_time_dimension(false);

  e_ti_.add_dimension("phi", phi_data_.length(), 
                      phi_data_.length(), 0);
  e_ti_.add_dimension("theta", theta_data_.length(), 
                      theta_data_.length(), 0);
  e_ti_.add_dimension("freq", frequencies_.length(),
                      frequencies_.length(), 0);
  e_ti_.set_name(base_name_ + "Et_imag");
  e_ti_.set_ptr(e_theta_imag_);
  e_ti_.has_time_dimension(false);

  e_pr_.add_dimension("phi", phi_data_.length(), 
                      phi_data_.length(), 0);
  e_pr_.add_dimension("theta", theta_data_.length(), 
                      theta_data_.length(), 0);
  e_pr_.add_dimension("freq", frequencies_.length(),
                      frequencies_.length(), 0);
  e_pr_.set_name(base_name_ + "Ep_real");
  e_pr_.set_ptr(e_phi_real_);
  e_pr_.has_time_dimension(false);

  e_pi_.add_dimension("phi", phi_data_.length(), 
                      phi_data_.length(), 0);
  e_pi_.add_dimension("theta", theta_data_.length(), 
                      theta_data_.length(), 0);
  e_pi_.add_dimension("freq", frequencies_.length(),
                      frequencies_.length(), 0);
  e_pi_.set_name(base_name_ + "Ep_imag");
  e_pi_.set_ptr(e_phi_imag_);
  e_pi_.has_time_dimension(false);

  // RCS
  rcs_var_.add_dimension("phi", phi_data_.length(), 
                      phi_data_.length(), 0);
  rcs_var_.add_dimension("theta", theta_data_.length(), 
                      theta_data_.length(), 0);
  rcs_var_.add_dimension("freq", frequencies_.length(),
                      frequencies_.length(), 0);
  rcs_var_.set_name(base_name_ + "rcs");
  rcs_var_.set_ptr(rcs_);
  rcs_var_.has_time_dimension(false);
  
  if (MPI_RANK == 0)
  {
    phi_.set_num(phi_data_.length());
    theta_.set_num(theta_data_.length());
    freqs_.set_num(frequencies_.length());

    E_theta_var_.set_num(sz);
    E_phi_var_.set_num(sz);

    e_tr_.set_num(dftsz);
    e_ti_.set_num(dftsz);
    e_pr_.set_num(dftsz);
    e_pi_.set_num(dftsz);
    rcs_var_.set_num(dftsz);
  }
  else
  {
    phi_.set_num(0);
    theta_.set_num(0);
    freqs_.set_num(0);

    E_theta_var_.set_num(0);
    E_phi_var_.set_num(0);
  }

  // TESTING ONLY!!!
  idx_tests();

}
  
void FarfieldResult::deinit()
{ 
  if (Wx_)
  {
    delete[] Wx_;
    Wx_ = 0;
  }

  if (Wy_)
  {
    delete[] Wy_;
    Wy_ = 0;
  }

  if (Wz_)
  {
    delete[] Wz_;
    Wz_ = 0;
  }

  if (Ux_)
  {
    delete[] Ux_;
    Ux_ = 0;
  }

  if (Uy_)
  {
    delete[] Uy_;
    Uy_ = 0;
  }

  if (Uz_)
  {
    delete[] Uz_;
    Uz_ = 0;
  }

  if (E_theta_)
  {
    delete[] E_theta_;
    E_theta_ = 0;
  }

  if (E_phi_)
  {
    delete[] E_phi_;
    E_phi_ = 0;
  }

  if (e_theta_real_)
  {
    delete[] e_theta_real_;
    delete[] e_theta_imag_;
    delete[] e_phi_real_;
    delete[] e_phi_imag_;
    delete[] rcs_;

    e_theta_real_ = 0;
    e_theta_imag_ = 0;
    e_phi_real_ = 0;
    e_phi_imag_ = 0;
    rcs_ = 0;
  }
}

ostream& FarfieldResult::to_string(ostream &os) const
{
  return os << "Near to Far field transformation using Luebbers' method.";
}

/**
 * Contains any data that is required by the FFAlg::alg() function
 */ 
class FFData
{
public:
  // W and U in farfield. These must be set to the correct pointers
  // for the current theta and phi for the observation point of interest. 
  // Thus, de-reference the pointer using the ff time step. 
  field_t *W_t1, *W_t2;
  field_t *U_t1, *U_t2;

  // An index into the above variables. 
  int idx;

  // Current time step
  unsigned int tstep;

  // Time delta
  delta_t dt;

  // Space deltas
  delta_t dx, dy, dz;

  // Centre of the grid
  grid_point grid_centre;

  // Shift along each axis to the location of the averaged value
  // Will be either 0.0 or 0.5, depending on the plane of intergration
  field_t xshift, yshift, zshift;

  // Distance from origin to observation sphere
  field_t obs_radius;

  // Cell size
  field_t cellsize;

  // # time to propagate between the two points which are
  // # farthest from each other on the Huygen's surface.
  field_t t_cross;

  unsigned int ff_tsteps;

  // Sign change, to make this match Jan's implementation. This is due
  // to the fact that we're approximating a time derivative across two
  // time steps, and the fact that the currents have various signs in
  // them.... not sure exactly where this comes from though...
  // This should be either -1 or +1. 
  field_t signchange;

  // Angle data
  Interval<field_t> *theta_data;
  Interval<field_t> *phi_data;

  FFData()
    : W_t1(0), W_t2(0), U_t1(0), U_t2(0),
      theta_data(0), phi_data(0)
  {}
};

/**
 * Compute the running sums for U and W at each time step and each
 * point of interest in the farfield.
 */
class FFAlg
{
public:
  static inline void alg(const int &x, const int &y, const int &z,
                         Fields_t &f, FFData &data)
  {
    field_t xt = (x - static_cast<int>(data.grid_centre.x)) * data.dx
      + (data.xshift * data.dx);
    field_t yt = (y - static_cast<int>(data.grid_centre.y)) * data.dy
      + (data.yshift * data.dy);
    field_t zt = (z - static_cast<int>(data.grid_centre.z)) * data.dz
      + (data.zshift * data.dz);

    // Distance to source point
    field_t r_prime = sqrt(xt*xt + yt*yt + zt*zt);

    // Angles to the far field observation point
    field_t phi, theta;
    int wuidx = 0;

    for (int phi_idx = 0; phi_idx < data.phi_data->length(); phi_idx++)
    {
      for (int theta_idx = 0; theta_idx < data.theta_data->length(); 
           theta_idx++)
      {
        phi = data.phi_data->get(phi_idx);
        theta = data.theta_data->get(theta_idx);
        wuidx = WU_index(phi_idx, theta_idx, 0, 
                         data.theta_data->length(), 
                         data.ff_tsteps);

        // double precision is required here because for floats,
        // floor(335.9999967972875) = 336!!!

        // Time is takes this source point to influence the farfield.
//     field_t tau = (data.obs_radius - (xt * sin(data.phi) * cos(data.theta)
//       + yt * sin(data.phi) * sin(data.theta) + zt * cos(data.phi))) / C;
        double tau = (data.tstep * data.dt) 
          - ((xt * sin(theta) * cos(phi)
              + yt * sin(theta) * sin(phi) 
              + zt * cos(theta))) / C
          + data.t_cross;

        // Farfield time step for this result. 
        double dtemp = tau / data.dt + 0.5;
        //int U_n = static_cast<int>(floor(tau / data.dt + 0.5));
        int U_n = static_cast<int>(floor(dtemp));
        dtemp = tau / data.dt;
        //int W_n = static_cast<int>(floor(tau / data.dt));
        int W_n = static_cast<int>(floor(dtemp));

        assert(U_n >= 0 && (U_n + 1) < data.ff_tsteps);
        assert(W_n >= 0 && (W_n + 1) < data.ff_tsteps);

        double U_a = (tau/data.dt + 0.5) - U_n;
        double W_a = (tau/data.dt) - W_n;

        assert(U_a >= 0.0 && U_a <= 1.0);
        assert(W_a >= 0.0 && W_a <= 1.0);

        field_t temp = data.signchange * data.cellsize; //  / (data.dt);

        data.U_t1[wuidx + U_n] += temp * (1 - U_a) * (f.et2_avg);
        data.U_t1[wuidx + U_n + 1] += temp * U_a * (f.et2_avg);

        data.U_t2[wuidx + U_n] -= temp * (1 - U_a) * (f.et1_avg);
        data.U_t2[wuidx + U_n + 1] -= temp * U_a * (f.et1_avg);
        
        data.W_t1[wuidx + W_n] -= (1 - W_a) * temp * (f.ht2_avg);
        data.W_t1[wuidx + W_n + 1] -= W_a * temp * (f.ht2_avg);

        data.W_t2[wuidx + W_n] += (1 - W_a) * temp * (f.ht1_avg);
        data.W_t2[wuidx + W_n + 1] += W_a * temp * (f.ht1_avg);
      }
    }

    data.idx++;
  }
};

map<string, Variable *> & FarfieldResult::get_result(const Grid &grid, 
                                                     unsigned int time_step)
{
  if (result_time(time_step))
  {
    // Set up the common stuff in the Data object
    FFData data;
    data.tstep = time_step;
    data.dt = grid.get_deltat();
    data.dx = grid.get_deltax();
    data.dy = grid.get_deltay();
    data.dz = grid.get_deltaz();
    data.grid_centre = grid.get_global_cell(grid.get_centre());
    data.obs_radius = r_;
    data.ff_tsteps = ff_tsteps_;
    data.t_cross = t_cross_;
    data.phi_data = &phi_data_;
    data.theta_data = &theta_data_;

    // Calculate W and U for each observation point of interest. 
    for (int face_idx = 0; face_idx < 6; face_idx++)
    {
      if (!use_face_[face_idx])
        continue;

      // This index is for the pointers storing the previous value
      // of the fields on the face.
      data.idx = 0;
      data.xshift = 0.5;
      data.yshift = 0.5;
      data.zshift = 0.5;

      switch (face_idx)
      {
      case FRONT:
      case BACK:
        data.W_t1 = Wy_;
        data.W_t2 = Wz_;
        data.U_t1 = Uy_;
        data.U_t2 = Uz_;

        data.cellsize = grid.get_deltay() * grid.get_deltaz();
        data.xshift = 0.0;
        break;
            
      case LEFT:
      case RIGHT:
        data.W_t1 = Wz_;
        data.W_t2 = Wx_;
        data.U_t1 = Uz_;
        data.U_t2 = Ux_;
        
        data.cellsize = grid.get_deltax() * grid.get_deltaz();
        data.yshift = 0.0;
        break;
        
      case TOP:
      case BOTTOM:
        data.W_t1 = Wx_;
        data.W_t2 = Wy_;
        data.U_t1 = Ux_;
        data.U_t2 = Uy_;
        
        data.cellsize = grid.get_deltax() * grid.get_deltay();
        data.zshift = 0.0;
        break;
      }            

      if (face_idx == BACK || face_idx == LEFT || face_idx == BOTTOM)
        data.signchange = -1.0;
      else
        data.signchange = 1.0;

      Face face = static_cast<Face>(face_idx);
      PlaneTiling<FFAlg, FFData>::loop(grid, *region_, 
                                       face, data);
    } // for face

  } // if (result_time...)

  return variables_;
}

map<string, Variable *> &FarfieldResult::get_pre_result(const Grid &grid)
{
  return pre_vars_;
}

map<string, Variable *> &FarfieldResult::get_post_result(const Grid &grid)
{
  if (MPI_RANK != 0)
    return post_vars_;
  
  // Calculate E_theta, E_phi, rcs
  //field_t temp = 1 / (4 * PI * r_ * C); // Non-normalized
  field_t temp = 1 / (4 * PI * C * grid.get_deltat()); // Normalized
  field_t W_t, W_p, U_t, U_p;

  int idx = 0;

  for (int phi_idx = 0; phi_idx < phi_data_.length(); phi_idx++)
  {
    for (int theta_idx = 0; theta_idx < theta_data_.length(); theta_idx++)
    {
      for (int fft_idx = 0; fft_idx < ff_tsteps_; fft_idx++)
      {
        field_t phi = phi_data_.get(phi_idx);
        field_t theta = theta_data_.get(theta_idx);

        W_t = temp * (Wx_[idx] * cos(theta) * cos(phi)
                      + Wy_[idx] * cos(theta) * sin(phi)
                      - Wz_[idx] * sin(theta));

        W_p = temp * (- Wx_[idx] * sin(phi) + Wy_[idx] * cos(phi));

        U_t = temp * (Ux_[idx] * cos(theta) * cos(phi)
                      + Uy_[idx] * cos(theta) * sin(phi)
                      - Uz_[idx] * sin(theta)); 

        U_p = temp * (- Ux_[idx] * sin(phi) + Uy_[idx] * cos(phi));

        E_theta_[idx] = - ETA * W_t - U_p;
        E_phi_[idx] = - ETA * W_p + U_t;

        int dft_idx = dft_index(phi_idx, theta_idx, 0);
        for (int f_idx = 0; f_idx < frequencies_.length(); f_idx++)
        {
          field_t e_cos = cos(-2 * PI * frequencies_.get(f_idx) 
                              * fft_idx * grid.get_deltat());
          field_t e_sin = sin(-2 * PI * frequencies_.get(f_idx) 
                              * fft_idx * grid.get_deltat());

          e_phi_real_[dft_idx] += E_phi_[idx] * e_cos;
          e_phi_imag_[dft_idx] += -1 * E_phi_[idx] * e_sin;

          e_theta_real_[dft_idx] += E_theta_[idx] * e_cos;
          e_theta_imag_[dft_idx] += -1 * E_theta_[idx] * e_sin;

          rcs_[dft_idx] += 
            (e_theta_imag_[dft_idx] * e_theta_imag_[dft_idx]
             + e_theta_real_[dft_idx] * e_theta_real_[dft_idx]
             + e_phi_imag_[dft_idx] * e_phi_imag_[dft_idx]
             + e_phi_real_[dft_idx] * e_phi_real_[dft_idx]) 
            / (2*ETA);

          dft_idx++;
        }

        idx++;
      } // end for ff_tsteps

    } // end for theta
  } // end for phi

  return post_vars_;
}

