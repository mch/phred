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

#include "FarfieldResult2.hh"
#include "../Constants.hh"
#include "../Globals.hh"
#include "../GridPlane.hh"

#include "../PlaneTiling.hh"

#include <fstream>
using namespace std;

void FarfieldResult2::export_dfts()
{
  ofstream jt1_of, jt2_of, mt1_of, mt2_of;

  jt1_of.open("jt1.txt", ofstream::out);
  jt2_of.open("jt2.txt", ofstream::out);
  mt1_of.open("mt1.txt", ofstream::out);
  mt2_of.open("mt2.txt", ofstream::out);

  for (int face_idx = 0; face_idx < 6; face_idx++)
  {
    if (!(*region_).has_face_data(static_cast<Face>(face_idx)))
      continue;

    region_t cells;

    cells.xmin = (*region_).xmin();
    cells.ymin = (*region_).ymin();
    cells.zmin = (*region_).zmin();
    cells.xmax = (*region_).xmax();
    cells.ymax = (*region_).ymax();
    cells.zmax = (*region_).zmax();

    switch (face_idx)
    {
    case FRONT:
      cells.xmin = cells.xmax - 1;
      break;

    case BACK:
      cells.xmax = cells.xmin + 1;
      break;

    case LEFT:
      cells.ymax = cells.ymin + 1;
      break;

    case RIGHT:
      cells.ymin = cells.ymax - 1;
      break;

    case TOP:
      cells.zmin = cells.zmax - 1;
      break;

    case BOTTOM:
      cells.zmax = cells.zmin + 1;
      break;
    }

    complex<field_t> *Jt1 = Jt1_data_[face_idx];
    complex<field_t> *Jt2 = Jt2_data_[face_idx];

    complex<field_t> *Mt1 = Mt1_data_[face_idx];
    complex<field_t> *Mt2 = Mt2_data_[face_idx];

    unsigned int idx = 0;

    for (unsigned int f_idx = 0; f_idx < frequencies_.length(); f_idx++)
    {
      for (unsigned int i = cells.xmin; i < cells.xmax; i++)
      {
        for (unsigned int j = cells.ymin; j < cells.ymax; j++)
        {
          for (unsigned int k = cells.zmin; k < cells.zmax; k++, idx++)
          {
            jt1_of << frequencies_.get(f_idx) << " " 
                   << i << " " << j << " " << k << " " 
                   << Jt1[idx].real() << " " << Jt1[idx].imag() << "\n";
            jt2_of << frequencies_.get(f_idx) << " " 
                   << i << " " << j << " " << k << " " 
                   << Jt2[idx].real() << " " << Jt2[idx].imag() << "\n";
            mt1_of << frequencies_.get(f_idx) << " " 
                   << i << " " << j << " " << k << " " 
                   << Mt1[idx].real() << " " << Mt1[idx].imag() << "\n";
            mt2_of << frequencies_.get(f_idx) << " " 
                   << i << " " << j << " " << k << " " 
                   << Mt2[idx].real() << " " << Mt2[idx].imag() << "\n";
          }
        }
      }
    }
  }
  
  jt1_of.close();
  jt2_of.close();
  mt1_of.close();
  mt2_of.close();
}

FarfieldResult2::FarfieldResult2()
  : r_(100), e_theta_data_(0), e_phi_data_(0), h_theta_data_(0),
    h_phi_data_(0), rcs_data_(0)
{
#ifdef HAVE_COMPLEX
  for (int i = 0; i < 6; i++)
  {
    Jt1_data_[i] = 0;
    Jt2_data_[i] = 0;
    Mt1_data_[i] = 0;
    Mt2_data_[i] = 0;

  }
#endif
}

FarfieldResult2::~FarfieldResult2()
{
  deinit();
}

void FarfieldResult2::set_radius(field_t r)
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

void FarfieldResult2::set_theta_degrees(field_t theta_start, 
                                        field_t theta_stop, 
                                        unsigned int num_theta)
{
  set_theta(theta_start * (PI / 180), theta_stop * (PI / 180),
            num_theta);
}

void FarfieldResult2::set_theta(field_t theta_start, field_t theta_stop, 
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

void FarfieldResult2::set_phi_degrees(field_t phi_start, field_t phi_stop, 
                                      unsigned int num_phi)
{
  set_phi(phi_start * (PI / 180), phi_stop * (PI / 180),
          num_phi);
}

void FarfieldResult2::set_phi(field_t phi_start, field_t phi_stop, 
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

void FarfieldResult2::init(const Grid &grid)
{
  if (box_.get())
  {
    region_ = grid.get_local_region(*(box_.get()));
  } else {
    throw ResultException("SurfaceCurrentResult has no surface defined!");
  }

  shared_ptr<Block> gregion = grid.get_global_region(*(box_.get()));

  // Setup the arrays for temporary and RCS data
  unsigned int sz = frequencies_.length() * theta_data_.length()
    * phi_data_.length();

  if (MPI_RANK == 0)
  {
    rcs_data_ = new field_t[sz];

#ifdef HAVE_COMPLEX
    e_theta_data_ = new complex<field_t>[sz];
    e_phi_data_ = new complex<field_t>[sz];
    h_theta_data_ = new complex<field_t>[sz];
    h_phi_data_ = new complex<field_t>[sz];
#endif
  }

  for (int i = 0; i < 6; i++)
  {
    unsigned int face_size = frequencies_.length();

    switch (i)
    {
    case FRONT:
    case BACK:
      face_size *= (*region_).ylen() * (*region_).zlen();
      break;

    case LEFT:
    case RIGHT:
      face_size *= (*region_).xlen() * (*region_).zlen();
      break;

    case TOP:
    case BOTTOM:
      face_size *= (*region_).xlen() * (*region_).ylen();
      break;
    }

    if ((*region_).has_face_data(static_cast<Face>(i)))
    {
      Jt1_data_[i] = new complex<field_t>[face_size];
      Jt2_data_[i] = new complex<field_t>[face_size];
      Mt1_data_[i] = new complex<field_t>[face_size];
      Mt2_data_[i] = new complex<field_t>[face_size];
    }
  }

  rcs_.reset();
  freqs_.reset();
  theta_.reset();
  phi_.reset();

//   post_vars_["e_theta"] = &E_theta_;
//   post_vars_["e_phi"] = &E_phi_;
//   post_vars_["h_theta"] = &H_theta_;
//   post_vars_["h_phi"] = &H_phi_;
  post_vars_["rcs"] = &rcs_;

  pre_vars_["freqs"] = &freqs_;
  pre_vars_["theta"] = &theta_;
  pre_vars_["phi"] = &phi_;

  rcs_.add_dimension("theta", theta_data_.length(), theta_data_.length(), 0);
  rcs_.add_dimension("phi", phi_data_.length(), phi_data_.length(), 0);
  rcs_.add_dimension("freqs", frequencies_.length(), frequencies_.length(), 0);
  rcs_.set_name(base_name_ + "rcs");
  rcs_.set_ptr(rcs_data_);
  rcs_.has_time_dimension(false);

//   E_theta_.add_dimension("theta", theta_data_.length(), 
//                          theta_data_.length(), 0);
//   E_theta_.add_dimension("phi", phi_data_.length(), 
//                          phi_data_.length(), 0);
//   E_theta_.add_dimension("freqs", frequencies_.length(), 
//                          frequencies_.length(), 0);
//   E_theta_.set_name(base_name_ + "e_theta");
//   E_theta_.set_ptr(e_theta_data_);
//   E_theta_.has_time_dimension(false);

//   E_phi_.add_dimension("theta", theta_data_.length(), 
//                          theta_data_.length(), 0);
//   E_phi_.add_dimension("phi", phi_data_.length(), 
//                          phi_data_.length(), 0);
//   E_phi_.add_dimension("freqs", frequencies_.length(), 
//                          frequencies_.length(), 0);
//   E_phi_.set_name(base_name_ + "e_phi");
//   E_phi_.set_ptr(e_phi_data_);
//   E_phi_.has_time_dimension(false);

//   H_theta_.add_dimension("theta", theta_data_.length(), 
//                          theta_data_.length(), 0);
//   H_theta_.add_dimension("phi", phi_data_.length(), 
//                          phi_data_.length(), 0);
//   H_theta_.add_dimension("freqs", frequencies_.length(), 
//                          frequencies_.length(), 0);
//   H_theta_.set_name(base_name_ + "h_theta");
//   H_theta_.set_ptr(h_theta_data_);
//   H_theta_.has_time_dimension(false);

//   H_phi_.add_dimension("theta", theta_data_.length(), 
//                          theta_data_.length(), 0);
//   H_phi_.add_dimension("phi", phi_data_.length(), 
//                          phi_data_.length(), 0);
//   H_phi_.add_dimension("freqs", frequencies_.length(), 
//                          frequencies_.length(), 0);
//   H_phi_.set_name(base_name_ + "h_phi");
//   H_phi_.set_ptr(h_phi_data_);
//   H_phi_.has_time_dimension(false);

  freqs_.add_dimension("freqs", frequencies_.length(), frequencies_.length(), 0);
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

  if (MPI_RANK == 0)
  {
    phi_.set_num(phi_data_.length());
    theta_.set_num(theta_data_.length());
    freqs_.set_num(frequencies_.length());
    //   H_phi_.set_num(sz);
    //   H_theta_.set_num(sz);
    //   E_phi_.set_num(sz);
    //   E_theta_.set_num(sz);
    rcs_.set_num(sz);
  }
  else
  {
    phi_.set_num(0);
    theta_.set_num(0);
    freqs_.set_num(0);
    //   H_phi_.set_num(0);
    //   H_theta_.set_num(0);
    //   E_phi_.set_num(0);
    //   E_theta_.set_num(0);
    rcs_.set_num(0);
  }
}

void FarfieldResult2::deinit()
{
  for (int i = 0; i < 6; i++)
  {
    if (Jt1_data_[i])
    {
      delete[] Jt1_data_[i];
      Jt1_data_[i] = 0;
    }

    if (Jt2_data_[i])
    {
      delete[] Jt2_data_[i];
      Jt2_data_[i] = 0;
    }

    if (Mt1_data_[i])
    {
      delete[] Mt1_data_[i];
      Mt1_data_[i] = 0;
    }

    if (Mt2_data_[i])
    {
      delete[] Mt2_data_[i];
      Mt2_data_[i] = 0;
    }
  }

  if (rcs_data_)
  {
    delete[] rcs_data_;
    rcs_data_ = 0;
  }

#ifdef HAVE_COMPLEX
  if (e_theta_data_)
  {
    delete[] e_theta_data_;
    e_theta_data_ = 0;
  }

  if (e_phi_data_)
  {
    delete[] e_phi_data_;
    e_phi_data_ = 0;
  }

  if (h_theta_data_)
  {
    delete[] h_theta_data_;
    h_theta_data_ = 0;
  }

  if (h_phi_data_)
  {
    delete[] h_phi_data_;
    h_phi_data_ = 0;
  }

#endif
}

map<string, Variable *> &
FarfieldResult2::get_pre_result(const Grid &grid)
{
  return pre_vars_;
}

map<string, Variable *> &
FarfieldResult2::get_post_result(const Grid &grid)
{
  // Calculate the far field!
  complex<field_t> eta(ETA,0);

  // TEMPORARY: checking only. Output the DFT's of the currents 
  export_dfts();  

  unsigned int idx = 0;
  for (int theta_idx = 0; theta_idx < theta_data_.length(); theta_idx++)
  {
    for (int phi_idx = 0; phi_idx < phi_data_.length(); phi_idx++)
    {
      for (int freq_idx = 0; freq_idx < frequencies_.length(); 
           freq_idx++, idx++)
      {
        // Calculate N and L by integrating J and M respectivly over
        // the surface of the box.
        vecp_t p;
        field_t k = 2 * PI * frequencies_.get(freq_idx) * sqrt(MU_0 * EPS_0);
        calc_potentials(p, theta_data_.get(theta_idx),
                        phi_data_.get(phi_idx),
                        freq_idx, k, grid);

        // Each rank will do this on it's own data... must reduce the result
        // to rank 0. 

        // Calculate E_phi, E_theta, H_phi, H_theta from L, N. 
        complex<field_t> temp(0,-1 * r_ * k);
        temp = complex<field_t>(0,1) * k * exp(temp);
        temp = temp / complex<field_t>(4 * PI * r_);

        e_theta_data_[idx] = complex<field_t>(-1,0) * temp 
          * (p.L_phi + eta * p.N_theta);

        e_phi_data_[idx] = temp * (p.L_theta - eta * p.N_phi);

        h_theta_data_[idx] = temp 
          * (p.N_phi - p.L_theta / eta);

        h_phi_data_[idx] =  complex<field_t>(-1,0) * temp 
          * (p.N_theta - p.L_phi / eta);

        // Calculate rcs
        field_t temp3 = (k * k) / (32 * (PI * PI) * ETA * r_ * r_);
        complex<field_t> temp1, temp2;
        temp1 = p.L_phi + eta * p.N_theta;
        temp2 = p.L_theta - eta * p.N_phi;

        field_t abs1 = abs(temp1);
        field_t abs2 = abs(temp2);
        rcs_data_[idx] = temp3 * (abs1 * abs1 + abs2 * abs2);

      } // end for freq
    } // end for phi
  } // end for theta

  return post_vars_;
}

void FarfieldResult2::calc_potentials(vecp_t &p, const field_t &theta, 
                                      const field_t &phi, 
                                      const unsigned int &f_idx, 
                                      const field_t &k, 
                                      const Grid &grid)
{
  grid_point grid_centre = grid.get_global_cell(grid.get_centre());

  field_t dx = grid.get_deltax();
  field_t dy = grid.get_deltay();
  field_t dz = grid.get_deltaz();

  for (int face_idx = 0; face_idx < 6; face_idx++)
  {
    if (!(*region_).has_face_data(static_cast<Face>(face_idx)))
      continue;

    region_t cells;

    complex<field_t> *Jt1 = Jt1_data_[face_idx];
    complex<field_t> *Jt2 = Jt2_data_[face_idx];

    complex<field_t> *Mt1 = Mt1_data_[face_idx];
    complex<field_t> *Mt2 = Mt2_data_[face_idx];

    cells.xmin = (*region_).xmin();
    cells.ymin = (*region_).ymin();
    cells.zmin = (*region_).zmin();
    cells.xmax = (*region_).xmax();
    cells.ymax = (*region_).ymax();
    cells.zmax = (*region_).zmax();

    unsigned int index = 0;

    switch (face_idx)
    {
    case FRONT:
    case BACK:
      if (face_idx == BACK)
        cells.xmax = cells.xmin + 1;
      else
        cells.xmin = cells.xmax - 1;

      index = f_idx * (*region_).ylen() * (*region_).zlen();

      // Jt1 == Jy, Jt2 == Jz. Jx == 0
      for (int idx = cells.xmin; idx < cells.xmax; idx++)
      {
        for (int jdx = cells.ymin; jdx < cells.ymax; jdx++)
        {
          for (int kdx = cells.zmin; kdx < cells.zmax; kdx++, index++)
          {
            field_t xt = (idx - static_cast<int>(grid_centre.x)) * dx;
            field_t yt = (jdx - static_cast<int>(grid_centre.y)) * dy;
            field_t zt = (kdx - static_cast<int>(grid_centre.z)) * dz;

            // Distance to source point
            field_t r_prime = sqrt(xt*xt + yt*yt + zt*zt);

            // Angles to source point
            field_t sp_phi = acos(zt / r_prime);
            field_t sp_theta = acos(xt / (r_prime * sin(sp_phi)));

            if (isnan(sp_theta)) // can happen if zt == r_prime
            {
              cout << "FB: sp_theta is NaN for xt = " << xt
                   << ", yt " << yt << ", zt = " << zt << ", r_prime = "
                   << r_prime << ", and sp_phi = " << sp_phi << endl;
              sp_theta = 0;
            }

            complex<field_t> cos_t_cos_p(cos(sp_theta) * cos(sp_phi), 0);
            complex<field_t> cos_t_sin_p(cos(sp_theta) * sin(sp_phi), 0);
            complex<field_t> sin_t(sin(sp_theta), 0);
            complex<field_t> sin_p(sin(sp_phi), 0);
            complex<field_t> cos_p(cos(sp_phi), 0);


            // Exponential phase term, r' cos \psi = \vec{r}' \cdot \hat{r}
            // CHECK THIS!
             field_t exp_phase = xt * sin(theta) * cos(phi) 
               + yt * sin(theta) * sin(phi)
               + zt * cos(theta);
//             field_t exp_phase = xt * sin(phi) * cos(theta)
//               + yt * sin(phi) * sin(theta) 
//               + zt * cos(phi);

            complex<field_t> temp(0, k * exp_phase);
            temp = exp(temp) * complex<field_t>(dy * dz, 0);

            p.N_theta += (Jt1[index] * cos_t_sin_p 
                          - Jt2[index] * sin_t) * temp;

            p.N_phi += (Jt1[index] * cos_p) * temp;
            
            p.L_theta += (Mt1[index] * cos_t_sin_p 
                          - Mt2[index] * sin_t) * temp;
            
            p.L_phi += (Mt1[index] * cos_p) * temp;
          }
        }
      }
      
      break;

    case LEFT:
    case RIGHT:
      if (face_idx == LEFT)
        cells.ymax = cells.ymin + 1;
      else
        cells.ymin = cells.ymax - 1;
        
      index = f_idx * (*region_).xlen() * (*region_).zlen();

      // Jt1 == Jz, Jt2 == Jx. Jy == 0
      for (int idx = cells.xmin; idx < cells.xmax; idx++)
      {
        for (int jdx = cells.ymin; jdx < cells.ymax; jdx++)
        {
          for (int kdx = cells.zmin; kdx < cells.zmax; kdx++, index++)
          {
            field_t xt = (idx - static_cast<int>(grid_centre.x)) * dx;
            field_t yt = (jdx - static_cast<int>(grid_centre.y)) * dy;
            field_t zt = (kdx - static_cast<int>(grid_centre.z)) * dz;

            // Distance to source point
            field_t r_prime = sqrt(xt*xt + yt*yt + zt*zt);

            // Angles to source point
            field_t sp_phi = acos(zt / r_prime);
            field_t sp_theta = acos(xt / (r_prime * sin(sp_phi)));

            if (isnan(sp_theta)) // can happen if zt == r_prime
            {
              cout << "LR: sp_theta is NaN for xt = " << xt
                   << ", yt " << yt << ", zt = " << zt << ", r_prime = "
                   << r_prime << ", and sp_phi = " << sp_phi << endl;
              sp_theta = 0;
            }

            complex<field_t> cos_t_cos_p(cos(sp_theta) * cos(sp_phi), 0);
            complex<field_t> cos_t_sin_p(cos(sp_theta) * sin(sp_phi), 0);
            complex<field_t> sin_t(sin(sp_theta), 0);
            complex<field_t> sin_p(sin(sp_phi), 0);
            complex<field_t> cos_p(cos(sp_phi), 0);

            // Exponential phase term, r' cos \psi = \vec{r}' \cdot \hat{r}
            // CHECK THIS!
            field_t exp_phase = xt * sin(theta) * cos(phi) 
              + yt * sin(theta) * sin(phi)
              + zt * cos(theta);

            complex<field_t> temp(0, k * exp_phase);

            temp = exp(temp) * complex<field_t>(dx * dz, 0);

            p.N_theta += (Jt2[index] * cos_t_cos_p
                          - Jt1[index] * sin_t) * temp;

            p.N_phi += (complex<field_t>(-1,0) * Jt2[index] * sin_p) * temp;
            
            p.L_theta += (Mt2[index] * cos_t_cos_p 
                          - Mt1[index] * sin_t) * temp;
            
            p.L_phi += (complex<field_t>(-1,0) * Mt2[index] * sin_p) * temp;
          }
        }
      }

      break;

    case TOP:
    case BOTTOM:
      if (face_idx == TOP)
        cells.zmin = cells.zmax - 1;
      else
        cells.zmax = cells.zmin + 1;

      index = f_idx * (*region_).ylen() * (*region_).xlen();

      // Jt1 == Jx, Jt2 == Jy. Jz == 0
      for (int idx = cells.xmin; idx < cells.xmax; idx++)
      {
        for (int jdx = cells.ymin; jdx < cells.ymax; jdx++)
        {
          for (int kdx = cells.zmin; kdx < cells.zmax; kdx++, index++)
          {
            field_t xt = (idx - static_cast<int>(grid_centre.x)) * dx;
            field_t yt = (jdx - static_cast<int>(grid_centre.y)) * dy;
            field_t zt = (kdx - static_cast<int>(grid_centre.z)) * dz;

            // Distance to source point
            field_t r_prime = sqrt(xt*xt + yt*yt + zt*zt);

            // Angles to source point
            field_t sp_phi = acos(zt / r_prime);
            field_t sp_theta = acos(xt / (r_prime * sin(sp_phi)));
            
            if (isnan(sp_theta)) // can happen if xt == 0 or yt == 0
            {
              sp_theta = 0;
            }

            complex<field_t> cos_t_cos_p(cos(sp_theta) * cos(sp_phi), 0);
            complex<field_t> cos_t_sin_p(cos(sp_theta) * sin(sp_phi), 0);
            complex<field_t> sin_t(sin(sp_theta), 0);
            complex<field_t> sin_p(sin(sp_phi), 0);
            complex<field_t> cos_p(cos(sp_phi), 0);

            // Exponential phase term, r' cos \psi = \vec{r}' \cdot \hat{r}
            // CHECK THIS!
            field_t exp_phase = xt * sin(theta) * cos(phi) 
              + yt * sin(theta) * sin(phi)
              + zt * cos(theta);

            complex<field_t> temp(0, k * exp_phase);

            temp = exp(temp) * complex<field_t>(dy * dx, 0);

            p.N_theta += (Jt1[index] * cos_t_cos_p
                          + Jt2[index] * cos_t_sin_p) * temp;

            p.N_phi += (complex<field_t>(-1,0) * Jt1[index] * sin_p
                        + Jt2[index] * cos_p) * temp;
            
            p.L_theta += (Mt1[index] * cos_t_cos_p 
                          + Mt2[index] * cos_t_sin_p) * temp;
            
            p.L_phi += (complex<field_t>(-1,0) * Mt1[index] * sin_p
                        + Mt2[index] * cos_p) * temp;

            // Ensure that there are no NaN's:
            bool have_nan = false;
            if (p.N_theta != p.N_theta)
            {
              cout << "p.N_theta is NaN! ";
              have_nan = true;
            }
            if (p.N_phi != p.N_phi)
            {
              cout << "p.N_phi is NaN! ";
              have_nan = true;
            }

            if (p.L_theta != p.L_theta)
            {
              cout << "p.N_theta is NaN! ";
              have_nan = true;
            }
            if (p.L_phi != p.L_phi)
            {
              cout << "p.L_phi is NaN! ";
              have_nan = true;
            }
            
            if (have_nan)
            {
              cout << ", temp = " << temp 
                   << ", index = " << index << ", Jt1[index] = "
                   << Jt1[index] << ", Jt2[index] = " << Jt2[index]
                   << ", Mt1[index] = " << Mt1[index] 
                   << ", Mt2[index] = " << Mt2[index] 
                   << ", idx = " << idx << ", jdx = " << jdx 
                   << ", kdx = " << kdx << ", theta = " << theta
                   << ", phi = " << phi << ", sp_theta = " 
                   << sp_theta << ", sp_phi = " << sp_phi << endl;
            }
            
          }
        }
      }


      break;
    }
    
  } // end for (face_idx...)
    
}

map<string, Variable *> &
FarfieldResult2::get_result(const Grid &grid, 
                            unsigned int time_step)
{
  // Just calculate data but return nothing until the end.
  for (int face_idx = 0; face_idx < 6; face_idx++)
  {
    if (!(*region_).has_face_data(static_cast<Face>(face_idx)))
      continue;

    region_t cells;

    cells.xmin = (*region_).xmin();
    cells.ymin = (*region_).ymin();
    cells.zmin = (*region_).zmin();
    cells.xmax = (*region_).xmax();
    cells.ymax = (*region_).ymax();
    cells.zmax = (*region_).zmax();

    switch (face_idx)
    {
    case FRONT:
      cells.xmin = cells.xmax - 1;
      calc_currents<YZPlane>(grid, time_step, cells,
                             face_idx);
      break;

    case BACK:
      cells.xmax = cells.xmin + 1;
      calc_currents<YZPlane>(grid, time_step, cells,
                             face_idx);
      break;

    case LEFT:
      cells.ymax = cells.ymin + 1;
      calc_currents<XZPlane>(grid, time_step, cells,
                             face_idx);
      break;

    case RIGHT:
      cells.ymin = cells.ymax - 1;
      calc_currents<XZPlane>(grid, time_step, cells,
                             face_idx);
      break;

    case TOP:
      cells.zmin = cells.zmax - 1;
      calc_currents<XYPlane>(grid, time_step, cells,
                             face_idx);
      break;

    case BOTTOM:
      cells.zmax = cells.zmin + 1;
      calc_currents<XYPlane>(grid, time_step, cells,
                             face_idx);
      break;
    }
  } // end for (face_idx...)

  return variables_;
}

template<class T>
void FarfieldResult2::calc_currents(const Grid &grid, 
                                    unsigned int time_step, 
                                    region_t &cells,
                                    int face_idx)
{
  field_t e_t1, e_t2, h_t1, h_t2;

  T p(const_cast<Grid &>(grid)); // EVIL

  complex<field_t> *Jt1 = Jt1_data_[face_idx];
  complex<field_t> *Jt2 = Jt2_data_[face_idx];

  complex<field_t> *Mt1 = Mt1_data_[face_idx];
  complex<field_t> *Mt2 = Mt2_data_[face_idx];

  unsigned int idx = 0;

  delta_t dt = grid.get_deltat();
  delta_t e_time = dt * time_step;
  delta_t h_time = dt * (static_cast<delta_t>(time_step) - 0.5);

  field_t e_cos_temp, e_sin_temp;
  field_t h_cos_temp, h_sin_temp;

  for (unsigned int f_idx = 0; f_idx < frequencies_.length(); f_idx++)
  {
    e_cos_temp = cos(-2 * PI * frequencies_.get(f_idx) * e_time);
    e_sin_temp = sin(-2 * PI * frequencies_.get(f_idx) * e_time);

    h_cos_temp = cos(-2 * PI * frequencies_.get(f_idx) * h_time);
    h_sin_temp = sin(-2 * PI * frequencies_.get(f_idx) * h_time);

    for (unsigned int i = cells.xmin; i < cells.xmax; i++)
    {
      for (unsigned int j = cells.ymin; j < cells.ymax; j++)
      {
        for (unsigned int k = cells.zmin; k < cells.zmax; k++, idx++)
        {
          e_t1 = p.get_avg_e_t1(i, j, k);
          e_t2 = p.get_avg_e_t2(i, j, k);
          h_t1 = p.get_avg_h_t1(i, j, k);
          h_t2 = p.get_avg_h_t2(i, j, k);

          // Slow due to temporaries? Optimized out? (On Kai C++,
          // apparently yes, optimized out)
          Jt1[idx] += complex<field_t>((-1) * h_t2 * h_cos_temp, 
                                       h_t2 * h_sin_temp);

          Jt2[idx] += complex<field_t>(h_t1 * h_cos_temp,
                                       (-1) * h_t1 * h_sin_temp);

          Mt1[idx] += complex<field_t>(e_t2 * e_cos_temp,
                                       (-1) * e_t2 * e_sin_temp);

          Mt2[idx] += complex<field_t>((-1) * e_t1 * e_cos_temp,
                                       e_t1 * e_sin_temp);

          // Above eqns verified, 2004-11-25 23:17 mch

          h_t2++; h_t1++; e_t1++; e_t2++;
        }
      }
    }
  } // end for f_idx
}

ostream& FarfieldResult2::to_string(ostream &os) const
{
#ifdef HAVE_COMPLEX
  return os << "FarfieldResult2; computing farfield...";
#else
  return os << "FarfieldResult2; complex types are not available. "
    "NOT CALCULATING FARFIELD.";
#endif
}

/**
 * A data object for use in the templated potential calculation
 * algorithms.
 */
class PotentialData
{
public:
  int idx;
    
  // Angles from x and z axis to observation point
  field_t theta, phi;

  // Frequency index and frequency in Hz
  int freq_idx;
  field_t freq;

  // Freespace wave number at the given frequency
  field_t k;

  // Origin of the grid; the point where the position vectors for
  // the observation and source points origionate. 
  point origin;
  grid_point grid_centre;

  // Cell sizes, since it will be necessary to convert cell numbers
  // to real coordinates
  field_t dx, dy, dz;

  // Area of a face of one cell on the plane of interest.
  field_t cell_area;

  // A object holding the potentials
  vecp_t p;

  // Pointers to the electric and magnetic current data that has
  // been accumulated on this face
  complex<field_t> *Jt1;
  complex<field_t> *Jt2;

  complex<field_t> *Mt1;
  complex<field_t> *Mt2;
};

/**
 * This algorithm is applied to each face of a Huygens surface (a box
 * in this case) and calculates N_theta, N_phi, L_theta, and L_phi by
 * integrating the electric and magnetic currents over the surface. 
 *
 * This class is templated by PotentialFunc, which must be a class
 * which implements a calc_potentials function. The reason for this is
 * that different potentials are required on different faces. 
 */ 
template<class PotentialFunc>
class PotentialAlg
{
public:
  static inline void alg(const int &x, const int &y, const int &z,
                         Fields_t &f, PotentialData &data)
  {
    field_t xt = (x - static_cast<int>(data.grid_centre.x)) * data.dx;
    field_t yt = (y - static_cast<int>(data.grid_centre.y)) * data.dy;
    field_t zt = (z - static_cast<int>(data.grid_centre.z)) * data.dz;

    // Distance to source point
    field_t r_prime = sqrt(xt*xt + yt*yt + zt*zt);

    // Angles to source point
    field_t sp_phi = acos(zt / r_prime);
    field_t sp_theta = acos(xt / (r_prime * sin(sp_phi)));
    
    complex<field_t> cos_t_cos_p(cos(sp_theta) * cos(sp_phi), 0);
    complex<field_t> cos_t_sin_p(cos(sp_theta) * sin(sp_phi), 0);
    complex<field_t> sin_t(sin(sp_theta), 0);
    complex<field_t> sin_p(sin(sp_phi), 0);
    complex<field_t> cos_p(cos(sp_phi), 0);
    
    // Exponential phase term, r' cos \psi = \vec{r}' \cdot \hat{r}
    // CHECK THIS! If this is correct, factor the sin and cos values
    // back to the data object. 
    field_t exp_phase = xt * sin(data.theta) * cos(data.phi) 
      + yt * sin(data.theta) * sin(data.phi)
      + zt * cos(data.theta);

    complex<field_t> temp(0, data.k * exp_phase);
    temp = exp(temp) * data.cell_area;

    // These depend on the face, because Jt1 and Jt2 are different for
    // each, but we need to know which is which so that the right
    // sin/cos constant can be multiplied in. 
    PotentialFunc::calc_potentials(cos_t_cos_p, cos_t_sin_p,
                                   sin_t, sin_p, cos_p, temp,
                                   data);

    data.idx++;
  }
};

class YZPotentials
{
public:
  static inline void calc_potentials(complex<field_t> &cos_t_cos_p,
                                     complex<field_t> &cos_t_sin_p,
                                     complex<field_t> &sin_t,
                                     complex<field_t> &sin_p,
                                     complex<field_t> &cos_p,
                                     complex<field_t> temp,
                                     PotentialData &data)
  {
    data.p.N_theta += (data.Jt1[data.idx] * cos_t_sin_p 
                  - data.Jt2[data.idx] * sin_t) * temp;
    
    data.p.N_phi += (data.Jt1[data.idx] * cos_p) * temp;
    
    data.p.L_theta += (data.Mt1[data.idx] * cos_t_sin_p 
                  - data.Mt2[data.idx] * sin_t) * temp;
    
    data.p.L_phi += (data.Mt1[data.idx] * cos_p) * temp;
  }
};

class XZPotentials
{
public:
  static inline void calc_potentials(complex<field_t> &cos_t_cos_p,
                                     complex<field_t> &cos_t_sin_p,
                                     complex<field_t> &sin_t,
                                     complex<field_t> &sin_p,
                                     complex<field_t> &cos_p,
                                     complex<field_t> temp,
                                     PotentialData &data)
  {
    data.p.N_theta += (data.Jt2[data.idx] * cos_t_cos_p
                  - data.Jt1[data.idx] * sin_t) * temp;
    
    data.p.N_phi += (complex<field_t>(-1,0) * data.Jt2[data.idx] 
                     * sin_p) * temp;
    
    data.p.L_theta += (data.Mt2[data.idx] * cos_t_cos_p 
                  - data.Mt1[data.idx] * sin_t) * temp;
    
    data.p.L_phi += (complex<field_t>(-1,0) * data.Mt2[data.idx] 
                     * sin_p) * temp;
  }
};

class XYPotentials
{
public:
  static inline void calc_potentials(complex<field_t> &cos_t_cos_p,
                                     complex<field_t> &cos_t_sin_p,
                                     complex<field_t> &sin_t,
                                     complex<field_t> &sin_p,
                                     complex<field_t> &cos_p,
                                     complex<field_t> temp,
                                     PotentialData &data)
  {
    data.p.N_theta += (data.Jt1[data.idx] * cos_t_cos_p
                  + data.Jt2[data.idx] * cos_t_sin_p) * temp;
    
    data.p.N_phi += (complex<field_t>(-1,0) * data.Jt1[data.idx] * sin_p
                + data.Jt2[data.idx] * cos_p) * temp;
    
    data.p.L_theta += (data.Mt1[data.idx] * cos_t_cos_p 
                  + data.Mt2[data.idx] * cos_t_sin_p) * temp;
    
    data.p.L_phi += (complex<field_t>(-1,0) * data.Mt1[data.idx] * sin_p
                + data.Mt2[data.idx] * cos_p) * temp;
  }
};
