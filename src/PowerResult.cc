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
#include <string.h> // for memset

PowerResult::PowerResult()
  : power_(0), step_x_(0), step_y_(0), step_z_(0), 
    plane_(0), normal_(X_AXIS)
{
  var_.has_time_dimension(false);
  freq_var_.has_time_dimension(false);
  var_.set_name("power");
  freq_var_.set_name("freqs");

  variables_["power"] = &var_;
  variables_["freqs"] = &freq_var_;
}

PowerResult::PowerResult(field_t freq_start, field_t freq_stop, 
                         unsigned int num_freqs)
  : DFTResult(freq_start, freq_stop, num_freqs), 
    power_(0), step_x_(0), step_y_(0), step_z_(0), 
    plane_(0), normal_(X_AXIS)
{
  var_.has_time_dimension(false);
  freq_var_.has_time_dimension(false);
  var_.set_name("power");
  freq_var_.set_name("freqs");

  variables_["power"] = &var_;
  variables_["freqs"] = &freq_var_;
}

~PowerResult::PowerResult()
{
  if (plane_)
    delete plane_;

  if (freqs_)
    delete[] freqs_;

  if (power_)
    delete[] power_;
}

void PowerResult::init(const Grid &grid)
{
  /* Region must be in out local sub-domain */ 
  region_ = grid.global_to_local(region_);

  /* Region must be a plane; set up grid plane */
  if (region_.xmax - region_.xmin == 1)
  {
    normal_ = X_AXIS;
    plane_ = new YZPlane();
    step_x_ = 1;
    cell_area_ = grid.get_deltay() * grid.get_deltaz();

  } else if (region_.ymax - region_.ymin == 1) {
    normal_ = Y_AXIS;
    plane_ = new XZPlane();
    step_y_ = -1;
    cell_area_ = grid.get_deltax() * grid.get_deltaz();

  } else if (region_.zmax - region_.zmin == 1) {
    normal_ = Z_AXIS;
    plane_ = new XYPlane();
    step_z_ = -1;
    cell_area_ = grid.get_deltax() * grid.get_deltay();

  } else {
    throw ResultException("Region must be limited to a plane for a PowerResult.");
  }

  /* Region must not be right up against the side of the domain
     because we need to average the H field with the next cell. */
 
  /* Set up the frequencies */ 
  freqs_ = new field_t[num_freqs_];
  power_ = new field_t[num_freqs_];

  memset(power_, 0, sizeof(field_t) * num_freqs_);

  freq_space = (freq_stop_ - freq_start_) / num_freqs_;

  for (int idx = 0; idx < num_freqs_; idx++)
    freqs_[idx] = freq_start_ + idx * freq_space_;

  /* Set up output variables. */
  var_.reset();
  freq_var_.reset();

  var_.add_dimension("Frequency", num_freqs_);
  freq_var_.add_dimension("Frequency", num_freqs_);
  var_.set_ptr(static_cast<void *>(power_));
  freq_var_.set_ptr(static_cast<void *>(freqs_));
  var_.set_num(1);
  freq_var_.set_num(1);
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

  if (power_)
  {
    delete[] power_;
    power_ = 0;
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
  field_t p_real;

  delta_t dt = grid.get_deltat();
  delta_t time = dt * time_step;

  field_t cos_temp, sin_temp;

#ifdef DEBUG
  cerr << "Computing power through a surface which is " 
       << region_.xmax - region_.xmin << "x" << 
       << region_.xmax - region_.xmin << "x" << 
       << region_.xmax - region_.xmin << " in size." << endl;
  cerr << "Frequency range: " << freq_start_ << " to " 
       << freq_stop_ << ", spacing: " << freq_space_ << ", number: "
       << num_freqs_ << endl;
  cerr << "Cell area is " << cell_area_ << endl;
#endif 

  for (unsigned int i = 0; i <= num_freqs_; i++)
  {
    cos_temp = cos(2 * PI * freqs_[i] * time);
    sin_temp = sin(2 * PI * freqs_[i] * time);
    
    p_real = 0;

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
          et1_imag = et1 * sin_temp;
          
          et2_real = et2 * cos_temp;
          et2_imag = et2 * sin_temp;
          
          ht1_real = ht1 * cos_temp;
          ht1_imag = ht1 * sin_temp;
          
          ht2_real = ht2 * cos_temp;
          ht2_imag = ht2 * sin_temp;

          p_real += (et1_real * ht2_real + et1_imag * ht2_imag)
            - (et2_real * ht1_real + et2_imag * ht1_imag);
        }
      }
    }
   
    power_[i] += 0.5 * p_real * cell_area_;
  }

  return variables_;
}
