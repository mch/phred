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

#include "PowerResult.hh"
#include "Constants.hh"

#include <string.h> // for memset
#include <math.h>

PowerResult::PowerResult()
  : power_real_(0), power_imag_(0), step_x_(0), step_y_(0), step_z_(0), 
    plane_(0), normal_(X_AXIS)
{
  real_var_.has_time_dimension(false);
  imag_var_.has_time_dimension(false);
  freq_var_.has_time_dimension(false);
  power_var_.has_time_dimension(true);

  variables_["real_power"] = &real_var_;
  variables_["imag_power"] = &imag_var_;
  variables_["freqs"] = &freq_var_;
  variables_["time_power"] = &power_var_;
}

PowerResult::PowerResult(field_t freq_start, field_t freq_stop, 
                         unsigned int num_freqs)
  : DFTResult(freq_start, freq_stop, num_freqs), 
    power_real_(0), power_imag_(0), step_x_(0), step_y_(0), step_z_(0), 
    plane_(0), normal_(X_AXIS)
{
  real_var_.has_time_dimension(false);
  imag_var_.has_time_dimension(false);
  freq_var_.has_time_dimension(false);
  power_var_.has_time_dimension(true);

  variables_["real_power"] = &real_var_;
  variables_["imag_power"] = &imag_var_;
  variables_["freqs"] = &freq_var_;
  variables_["time_power"] = &power_var_;
}

PowerResult::~PowerResult()
{
  if (plane_)
    delete plane_;

  if (freqs_)
    delete[] freqs_;

  if (power_real_)
    delete[] power_real_;

  if (power_imag_)
    delete[] power_imag_;

  if (et1r_)
  {
    delete[] et1r_;
    delete[] et1i_;
    delete[] et2r_;
    delete[] et2i_;
    delete[] ht1r_;
    delete[] ht1i_;
    delete[] ht2r_;
    delete[] ht2i_;

    et1r_ = 0;
    et1i_ = 0;
    et2r_ = 0;
    et2i_ = 0;

    ht1r_ = 0;
    ht1i_ = 0;
    ht2r_ = 0;
    ht2i_ = 0;
  }
}

void PowerResult::init(const Grid &grid)
{
  /* Region must be in out local sub-domain */ 
  region_ = grid.global_to_local(region_);
  x_size_ = region_.xmax - region_.xmin;
  y_size_ = region_.ymax - region_.ymin;
  z_size_ = region_.zmax - region_.zmin;

  /* Region must be a plane; set up grid plane */
  if (region_.xmax - region_.xmin == 1)
  {
    normal_ = X_AXIS;
    plane_ = new YZPlane(const_cast<Grid&>(grid));
    step_x_ = 1;
    cell_area_ = grid.get_deltay() * grid.get_deltaz();

  } else if (region_.ymax - region_.ymin == 1) {
    normal_ = Y_AXIS;
    plane_ = new XZPlane(const_cast<Grid&>(grid));
    step_y_ = -1;
    cell_area_ = grid.get_deltax() * grid.get_deltaz();

  } else if (region_.zmax - region_.zmin == 1) {
    normal_ = Z_AXIS;
    plane_ = new XYPlane(const_cast<Grid&>(grid));
    step_z_ = -1;
    cell_area_ = grid.get_deltax() * grid.get_deltay();

  } else {
    throw ResultException("Region must be limited to a plane for a PowerResult.");
  }

  /* GRR, ARGH! */
  unsigned int sz = x_size_ * y_size_ * z_size_;

  et1r_ = new field_t[sz];
  et1i_ = new field_t[sz];
  ht1r_ = new field_t[sz];
  ht1i_ = new field_t[sz];

  et2r_ = new field_t[sz];
  et2i_ = new field_t[sz];
  ht2r_ = new field_t[sz];
  ht2i_ = new field_t[sz];

  if (!et1r_ || !et1i_ || !ht1r_ || !ht1i_
      || !et2r_ || !et2i_ || !ht2r_ || !ht2i_)
    throw MemoryException();

  memset(et1r_, 0, sizeof(field_t) * sz);
  memset(et1i_, 0, sizeof(field_t) * sz);
  memset(et2r_, 0, sizeof(field_t) * sz);
  memset(et2i_, 0, sizeof(field_t) * sz);

  memset(ht1r_, 0, sizeof(field_t) * sz);
  memset(ht1i_, 0, sizeof(field_t) * sz);
  memset(ht2r_, 0, sizeof(field_t) * sz);
  memset(ht2i_, 0, sizeof(field_t) * sz);

  /* Region must not be right up against the side of the domain
     because we need to average the H field with the next cell. */
 
  /* Set up the frequencies */ 
  freqs_ = new field_t[num_freqs_ + 1];
  power_real_ = new field_t[num_freqs_ + 1];
  power_imag_ = new field_t[num_freqs_ + 1];

  if (!freqs_ || !power_real_ || !power_imag_)
    throw MemoryException();

  memset(power_imag_, 0, sizeof(field_t) * (num_freqs_ + 1));
  memset(power_real_, 0, sizeof(field_t) * (num_freqs_ + 1));

  freq_space_ = (freq_stop_ - freq_start_) / num_freqs_;

  for (int idx = 0; idx <= num_freqs_; idx++)
    freqs_[idx] = freq_start_ + idx * freq_space_;

  /* Set up output variables. */
  real_var_.reset();
  imag_var_.reset();
  freq_var_.reset();
  power_var_.reset();

  real_var_.set_name(base_name_ + "_power_real");
  imag_var_.set_name(base_name_ + "_power_imag");
  freq_var_.set_name(base_name_ + "_freqs");
  power_var_.set_name(base_name_ + "_time_power");

  real_var_.add_dimension("Frequency", num_freqs_);
  imag_var_.add_dimension("Frequency", num_freqs_);
  freq_var_.add_dimension("Frequency", num_freqs_);
  power_var_.add_dimension("Power", 1, 0);

  imag_var_.set_ptr(power_imag_);
  real_var_.set_ptr(power_real_);  
  freq_var_.set_ptr(freqs_);
  power_var_.set_ptr(&time_power_);
  real_var_.set_num(num_freqs_);
  imag_var_.set_num(num_freqs_);
  freq_var_.set_num(num_freqs_);
  power_var_.set_num(1);
}

void PowerResult::deinit(const Grid &grid)
{
  if (plane_)
  {
    delete plane_;
    plane_ = 0;
  }

  if (freqs_)
  {
    delete[] freqs_;
    freqs_ = 0;
  }

  if (power_real_)
  {
    delete[] power_real_;
    power_real_ = 0;
  }

  if (power_imag_)
  {
    delete[] power_imag_;
    power_imag_ = 0;
  }

  if (et1r_)
  {
    delete[] et1r_;
    delete[] et1i_;
    delete[] et2r_;
    delete[] et2i_;
    delete[] ht1r_;
    delete[] ht1i_;
    delete[] ht2r_;
    delete[] ht2i_;

    et1r_ = 0;
    et1i_ = 0;
    et2r_ = 0;
    et2i_ = 0;

    ht1r_ = 0;
    ht1i_ = 0;
    ht2r_ = 0;
    ht2i_ = 0;
  }
}

map<string, Variable *> &PowerResult::get_result(const Grid &grid, 
                                                 unsigned int time_step)
{
  field_t et1, et2, ht1, ht2;
  field_t et1_real, et1_imag;
  field_t et2_real, et2_imag;
  field_t ht1_real, ht1_imag;
  field_t ht2_real, ht2_imag;
  field_t p_real2, p_imag2;

  delta_t dt = grid.get_deltat();
  delta_t time = dt * time_step;

  field_t cos_temp, sin_temp;

// #ifdef DEBUG
//   if (time_step == 10 || time_step == 15)
//   {
//     cerr << "Computing power through a surface which is "
//          << region_.xmax - region_.xmin << "x"
//          << region_.ymax - region_.ymin << "x"
//          << region_.zmax - region_.zmin << " in size." << endl;
//     cerr << "Frequency range: " << freq_start_ << " to " 
//          << freq_stop_ << ", spacing: " << freq_space_ << ", number: "
//          << num_freqs_ << endl;
//     cerr << "Cell area is " << cell_area_ << endl;
//   }
// #endif 

  // Rearrange this so that there is only one loop traversing the
  // data, and put the time and freq calculations inside?

  // Compute the instantaneous power through the surface at this
  // instant in time.
  unsigned int idx = 0;
  time_power_ = 0;

  for (unsigned int x = region_.xmin; x < region_.xmax; x++) 
  {
    for (unsigned int y = region_.ymin; y < region_.ymax; y++) 
    {
      for (unsigned int z = region_.zmin; z < region_.zmax; z++) 
      {
        et1 = plane_->get_e_t1(x, y, z);
        et2 = plane_->get_e_t2(x, y, z);
        
        ht1 = (plane_->get_h_t1(x, y, z) 
               + plane_->get_h_t1(x + step_x_, y + step_y_, 
                                  z + step_z_)) / 2;
        ht2 = (plane_->get_h_t2(x, y, z) 
               + plane_->get_h_t2(x + step_x_, y + step_y_, 
                                  z + step_z_)) / 2;

        time_power_ += et1 * ht2 - et2 * ht1;
      }
    }
  }

  // Compute the power throught the surface in the frequency domain
  for (unsigned int i = 0; i <= num_freqs_; i++)
  {
    cos_temp = cos(2 * PI * freqs_[i] * time);
    sin_temp = sin(2 * PI * freqs_[i] * time);

    p_real2 = 0;
    p_imag2 = 0;

    unsigned int idx = 0;

    for (unsigned int x = region_.xmin; x < region_.xmax; x++) 
    {
      for (unsigned int y = region_.ymin; y < region_.ymax; y++) 
      {
        for (unsigned int z = region_.zmin; z < region_.zmax; z++) 
        {
          et1 = plane_->get_e_t1(x, y, z);
          et2 = plane_->get_e_t2(x, y, z);

          ht1 = (plane_->get_h_t1(x, y, z) 
                 + plane_->get_h_t1(x + step_x_, y + step_y_, 
                                    z + step_z_)) / 2;
          ht2 = (plane_->get_h_t2(x, y, z) 
                 + plane_->get_h_t2(x + step_x_, y + step_y_, 
                                    z + step_z_)) / 2;
          
          et1_real = et1 * cos_temp;
          et1_imag = (-1) * et1 * sin_temp;
          
          et2_real = et2 * cos_temp;
          et2_imag = (-1) * et2 * sin_temp;
          
          ht1_real = ht1 * cos_temp;
          ht1_imag = (-1) * ht1 * sin_temp;
          
          ht2_real = ht2 * cos_temp;
          ht2_imag = (-1) * ht2 * sin_temp;
          
          et1r_[idx] += et1_real;
          et2r_[idx] += et2_real;
          et1i_[idx] += et1_imag;
          et2i_[idx] += et2_imag;

          ht1r_[idx] += ht1_real;
          ht2r_[idx] += ht2_real;
          ht1i_[idx] += ht1_imag;
          ht2i_[idx] += ht2_imag;

          p_real2 += (et1r_[idx] * ht2r_[idx] + et1i_[idx] * ht2i_[idx])
            - (et2r_[idx] * ht1r_[idx] + et2i_[idx] * ht1i_[idx]);

          p_imag2 += et1i_[idx] * ht2r_[idx] - ht2i_[idx] * et1r_[idx] 
            + ht1i_[idx] * et2r_[idx] - et2i_[idx] * ht1r_[idx];

          ++idx;
        }
      }
    }
   
    power_real_[i] = 0.5 * p_real2 * cell_area_;
    power_imag_[i] = 0.5 * p_imag2 * cell_area_;

  }

  return variables_;
}