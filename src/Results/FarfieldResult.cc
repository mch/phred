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
          << WU_index(phi_idx, theta_idx, 0) << "\n";
    }
  }

  wuf.close();

  ofstream tf;
  tf.open("tempidxtest.txt", ofstream::out);

  tf << "faceidx comp temp_idx(0,0,10,10)\n";

  for (int faceidx = 0; faceidx < 6; faceidx++)
  {
    for (int comp = 0; comp < 3; comp++)
    {
      tf << faceidx << "\t" << comp << "\t" 
         << temp_index(faceidx, comp, 0,0, 10, 10) << "\n";
    }
  }

  tf.close();
}

FarfieldResult::FarfieldResult()
  : r_(100),
    Wx_(0), Wy_(0), Wz_(0), Ux_(0), Uy_(0), Uz_(0),
    E_temp_(0), H_temp_(0), E_theta_(0), E_phi_(0)
{}

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

  t_cross_ = static_cast<unsigned int>
    (ceil(maxlength / (C * grid.get_deltat())));

  ff_tsteps_ = (time_stop_ - time_start_) + 2 * t_cross_;
  
#ifdef DEBUG
  cerr << "Farfield, num ff_tsteps: " << ff_tsteps_ << endl;
#endif

  // Allocate space for W and U
  unsigned int sz = ff_tsteps_ * theta_data_.length() * phi_data_.length();

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
  }

  // Allocate space to store past values of E and H on the faces so
  // that the time derivate can be properly approximated.
  unsigned int tempsz = 2 * region_->xlen() * region_->ylen()
    + 2 * region_->xlen() * region_->zlen()
    + 2 * region_->ylen() * region_->zlen();
  
  E_temp_ = new field_t[tempsz * 6 * 3];
  H_temp_ = new field_t[tempsz * 6 * 3];

  // Set up output variables
  freqs_.reset();
  theta_.reset();
  phi_.reset();
  E_theta_var_.reset();
  E_phi_var_.reset();

  pre_vars_["freqs"] = &freqs_;
  pre_vars_["theta"] = &theta_;
  pre_vars_["phi"] = &phi_;

  post_vars_["Et"] = &E_theta_var_;
  post_vars_["Ep"] = &E_phi_var_;

  freqs_.add_dimension("freqs", frequencies_.length(), 
                       frequencies_.length(), 0);
  freqs_.set_name(base_name_ + "freqs");
  freqs_.set_ptr(frequencies_.get_ptr());
  freqs_.has_time_dimension(false);

  theta_.add_dimension("theta", theta_data_.length(), 
                       theta_data_.length(), 0);
  theta_.set_name(base_name_ + "theta");
  theta_.set_ptr(theta_data_.get_ptr());
  theta_.has_time_dimension(false);

  phi_.add_dimension("phi", phi_data_.length(), phi_data_.length(), 0);
  phi_.set_name(base_name_ + "phi");
  phi_.set_ptr(phi_data_.get_ptr());
  phi_.has_time_dimension(false);

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

  if (MPI_RANK == 0)
  {
    phi_.set_num(phi_data_.length());
    theta_.set_num(theta_data_.length());
    freqs_.set_num(frequencies_.length());

    E_theta_var_.set_num(sz);
    E_phi_var_.set_num(sz);
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

  if (E_temp_)
  {
    delete[] E_temp_;
    E_temp_ = 0;
  }

  if (H_temp_)
  {
    delete[] H_temp_;
    H_temp_ = 0;
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

  // E and H temp so that we can properly calculate the derivate of
  // the components w.r.t. time. These must be dereferenced to the
  // correct face and field components.
  field_t *E_t1, *E_t2, *H_t1, *H_t2;

  // An index into the above variables. 
  int idx;

  // Current point in the far field for which W and U are being
  // evaluated.
  field_t theta, phi;

  // Current time step
  unsigned int tstep;

  // Time delta
  delta_t dt;

  // Space deltas
  delta_t dx, dy, dz;

  // Centre of the grid
  grid_point grid_centre;

  // Distance from origin to observation sphere
  field_t obs_radius;

  // Cell size
  field_t cellsize;

  // # timesteps to propagate between the two points which are
  // # farthest from each other on the Huygen's surface.
  unsigned int t_cross;

  // Debugging only
  unsigned int ff_tsteps_;

  FFData()
    : W_t1(0), W_t2(0), U_t1(0), U_t2(0), E_t1(0), E_t2(0), 
      H_t1(0), H_t2(0)
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
    field_t xt = (x - static_cast<int>(data.grid_centre.x)) * data.dx;
    field_t yt = (y - static_cast<int>(data.grid_centre.y)) * data.dy;
    field_t zt = (z - static_cast<int>(data.grid_centre.z)) * data.dz;

    // Distance to source point
    field_t r_prime = sqrt(xt*xt + yt*yt + zt*zt);

    // Time is takes this source point to influence the farfield.
//     field_t tau = (data.obs_radius - (xt * sin(data.phi) * cos(data.theta)
//       + yt * sin(data.phi) * sin(data.theta) + zt * cos(data.phi))) / C;
    field_t tau = ((xt * sin(data.phi) * cos(data.theta)
      + yt * sin(data.phi) * sin(data.theta) + zt * cos(data.phi))) 
      / (C * data.dt);

    tau = tau + data.t_cross;

    // Farfield time step for this result. 
    int U_n = static_cast<int>(floor(data.tstep + 0.5 + tau));
    int W_n = static_cast<int>(floor(data.tstep + tau));

    assert(U_n >= 0 && (U_n + 1) < data.ff_tsteps_);
    assert(W_n >= 0 && (W_n + 1) < data.ff_tsteps_);

    field_t U_a = (data.tstep + 0.5 + tau) - U_n;
    field_t W_a = (data.tstep + tau) - W_n;

    assert(U_a >= 0.0 && U_a <= 1.0);
    assert(W_a >= 0.0 && W_a <= 1.0);

    // The 4 * PI * r * C has been factored out and must be applied
    // when the spherical components are calculated.
    //field_t temp = data.cellsize / (4 * PI * data.obs_radius * C * data.dt);
    field_t temp = data.cellsize / (data.dt);

    data.U_t1[U_n] += temp * (1 - U_a) * (f.et2_avg - data.E_t2[data.idx]);
    data.U_t2[U_n] -= temp * (1 - U_a) * (f.et1_avg - data.E_t1[data.idx]);

    data.U_t1[U_n + 1] += temp * U_a * (f.et2_avg - data.E_t2[data.idx]);
    data.U_t2[U_n + 1] -= temp * U_a * (f.et1_avg - data.E_t1[data.idx]);

    data.E_t1[data.idx] = f.et1_avg;
    data.E_t2[data.idx] = f.et2_avg;

    data.W_t1[W_n] -= (1 - W_a) * temp * (f.ht2_avg - data.H_t2[data.idx]);
    data.W_t2[W_n] += (1 - W_a) * temp * (f.ht1_avg - data.H_t1[data.idx]);

    data.W_t1[W_n + 1] -= W_a * temp * (f.ht2_avg - data.H_t2[data.idx]);
    data.W_t2[W_n + 1] += W_a * temp * (f.ht1_avg - data.H_t1[data.idx]);

    data.H_t1[data.idx] = f.ht1_avg;
    data.H_t2[data.idx] = f.ht2_avg;

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
    data.ff_tsteps_ = ff_tsteps_;
    data.t_cross = t_cross_;

    // Calculate W and U for each observation point of interest. 
    for (int phi_idx = 0; phi_idx < phi_data_.length(); phi_idx++)
    {
      for (int theta_idx = 0; theta_idx < theta_data_.length(); theta_idx++)
      {
        data.phi = phi_data_.get(phi_idx);
        data.theta = theta_data_.get(theta_idx);
        
        for (int face_idx = 0; face_idx < 6; face_idx++)
        {
          data.idx = 0;
          
          switch (face_idx)
          {
          case FRONT:
          case BACK:
            data.W_t1 = &Wy_[WU_index(phi_idx, theta_idx, 0)];
            data.W_t2 = &Wz_[WU_index(phi_idx, theta_idx, 0)];
            data.U_t1 = &Uy_[WU_index(phi_idx, theta_idx, 0)];
            data.U_t2 = &Uz_[WU_index(phi_idx, theta_idx, 0)];

            data.E_t1 = &E_temp_[temp_index(face_idx, 1, 0, 0, 
                                            region_->ylen(), region_->zlen())];
            data.E_t2 = &E_temp_[temp_index(face_idx, 2, 0, 0, 
                                            region_->ylen(), region_->zlen())];

            data.H_t1 = &H_temp_[temp_index(face_idx, 1, 0, 0, 
                                            region_->ylen(), region_->zlen())];
            data.H_t2 = &H_temp_[temp_index(face_idx, 2, 0, 0, 
                                            region_->ylen(), region_->zlen())];

            data.cellsize = grid.get_deltay() * grid.get_deltaz();
            break;
            
          case LEFT:
          case RIGHT:
            data.W_t1 = &Wz_[WU_index(phi_idx, theta_idx, 0)];
            data.W_t2 = &Wx_[WU_index(phi_idx, theta_idx, 0)];
            data.U_t1 = &Uz_[WU_index(phi_idx, theta_idx, 0)];
            data.U_t2 = &Ux_[WU_index(phi_idx, theta_idx, 0)];

            data.E_t1 = &E_temp_[temp_index(face_idx, 2, 0, 0, 
                                          region_->ylen(), region_->zlen())];
            data.E_t2 = &E_temp_[temp_index(face_idx, 0, 0, 0, 
                                          region_->ylen(), region_->zlen())];
            data.H_t1 = &H_temp_[temp_index(face_idx, 2, 0, 0, 
                                            region_->ylen(), region_->zlen())];
            data.H_t2 = &H_temp_[temp_index(face_idx, 0, 0, 0, 
                                            region_->ylen(), region_->zlen())];

            data.cellsize = grid.get_deltax() * grid.get_deltaz();
            break;

          case TOP:
          case BOTTOM:
            data.W_t1 = &Wx_[WU_index(phi_idx, theta_idx, 0)];
            data.W_t2 = &Wy_[WU_index(phi_idx, theta_idx, 0)];
            data.U_t1 = &Ux_[WU_index(phi_idx, theta_idx, 0)];
            data.U_t2 = &Uy_[WU_index(phi_idx, theta_idx, 0)];

            data.E_t1 = &E_temp_[temp_index(face_idx, 0, 0, 0, 
                                          region_->ylen(), region_->zlen())];
            data.E_t2 = &E_temp_[temp_index(face_idx, 1, 0, 0, 
                                          region_->ylen(), region_->zlen())];

            data.H_t1 = &H_temp_[temp_index(face_idx, 0, 0, 0, 
                                            region_->ylen(), region_->zlen())];
            data.H_t2 = &H_temp_[temp_index(face_idx, 1, 0, 0, 
                                            region_->ylen(), region_->zlen())];

            data.cellsize = grid.get_deltax() * grid.get_deltay();
            break;
          }            

          PlaneTiling<FFAlg, FFData>::loop(grid, *region_, 
                                           static_cast<Face>(face_idx), 
                                           data);
        }
      } // end for theta
    } // end for phi
  }

  return variables_;
}

map<string, Variable *> &FarfieldResult::get_pre_result(const Grid &grid)
{
  return pre_vars_;
}

map<string, Variable *> &FarfieldResult::get_post_result(const Grid &grid)
{
  // Calculate E_theta, E_phi, rcs
  field_t temp = 4 * PI * r_ * C;
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
                      - Uz_[idx] * sin(phi)); 

        U_p = temp * (- Ux_[idx] * sin(phi) + Uy_[idx] * cos(phi));

        E_theta_[idx] = - ETA * W_t - U_p;
        E_phi_[idx] = - ETA * W_p + U_t;

        idx++;
      }
    } // end for theta
  } // end for phi

  return post_vars_;
}

