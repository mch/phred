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
#include "../Constants.hh"
#include "../Globals.hh"

#include <string.h> // for memset
#include <math.h>

PowerResult::PowerResult()
  : freqs_(0), power_real_(0), power_imag_(0), time_power_(0),
    has_data_(false), xmin_(0), ymin_(0), zmin_(0), 
    xmax_(0), ymax_(0), zmax_(0),
    step_x_(0), step_y_(0), step_z_(0), 
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
    freqs_(0), power_real_(0), power_imag_(0), time_power_(0),
    has_data_(false), step_x_(0), step_y_(0), step_z_(0), 
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
  const GridInfo &gi = grid.get_grid_info();

  /* Region must be in out local sub-domain */ 
  region_ = grid.get_local_region(*box_);

  x_size_ = (*region_).xmax() - (*region_).xmin();
  y_size_ = (*region_).ymax() - (*region_).ymin();
  z_size_ = (*region_).zmax() - (*region_).zmin();

  xmin_ = (*region_).xmin();
  ymin_ = (*region_).ymin();
  zmin_ = (*region_).zmin();

  xmax_ = (*region_).xmax();
  ymax_ = (*region_).ymax();
  zmax_ = (*region_).zmax();

  /* Region must be a plane; set up grid plane */
  if (face_ == FRONT || face_ == BACK)
  {
    normal_ = X_AXIS;
    plane_ = new YZPlane(const_cast<Grid&>(grid));
    step_x_ = 1;
    x_size_ = 1;
    cell_area_ = grid.get_deltay() * grid.get_deltaz();

    if (face_ == FRONT && (*region_).has_face_data(FRONT))
    {
      has_data_ = true;
      xmin_ = xmax_ - 1;
    }

    if (face_ == BACK && (*region_).has_face_data(BACK))
    {
      has_data_ = true;
      xmax_ = xmin_ + 1;
    }
  } 
  else if (face_ == LEFT || face_ == RIGHT) 
  {
    normal_ = Y_AXIS;
    plane_ = new XZPlane(const_cast<Grid&>(grid));
    step_y_ = -1;
    y_size_ = 1;
    cell_area_ = grid.get_deltax() * grid.get_deltaz();

    if (face_ == RIGHT && (*region_).has_face_data(RIGHT))
    {
      has_data_ = true;
      ymin_ = ymax_ - 1;
    }

    if (face_ == LEFT && (*region_).has_face_data(LEFT))
    {
      has_data_ = true;
      ymax_ = ymin_ + 1;
    }
  } 
  else if (face_ == TOP || face_ == BOTTOM)
  {
    normal_ = Z_AXIS;
    plane_ = new XYPlane(const_cast<Grid&>(grid));
    step_z_ = -1;
    z_size_ = 1;
    cell_area_ = grid.get_deltax() * grid.get_deltay();

    if (face_ == TOP && (*region_).has_face_data(TOP))
    {
      has_data_ = true;
      zmin_ = zmax_ - 1;
    }

    if (face_ == BOTTOM && (*region_).has_face_data(BOTTOM))
    {
      has_data_ = true;
      zmax_ = zmin_ + 1;
    }
  } else {
    throw ResultException("Region must be limited to a plane for a PowerResult.");
  }

  // Must always set up the variables, regardless of wether we
  // contribute any data, otherwise we may have problems with MPI in
  // the DataWriter.
  real_var_.reset();
  imag_var_.reset();
  freq_var_.reset();
  power_var_.reset();
  
  real_var_.set_name(base_name_ + "_power_real");
  imag_var_.set_name(base_name_ + "_power_imag");
  freq_var_.set_name(base_name_ + "_freqs");
  power_var_.set_name(base_name_ + "_time_power");
  
  real_var_.add_dimension("Frequency", num_freqs_, num_freqs_, 0);
  imag_var_.add_dimension("Frequency", num_freqs_, num_freqs_, 0);
  freq_var_.add_dimension("Frequency", num_freqs_, num_freqs_, 0);
  power_var_.add_dimension("Power", 1, 1, 0);
  
  if (has_data_)
  {
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
    imag_var_.set_ptr(power_imag_);
    real_var_.set_ptr(power_real_);  
    freq_var_.set_ptr(freqs_);
    power_var_.set_ptr(&time_power_);

    if (MPI_RANK == 0) // All data is collected to rank 0 by this
                       // result, rather than by the DataWriters.
    {
      power_var_.set_num(1);
      real_var_.set_num(num_freqs_);
      imag_var_.set_num(num_freqs_);
      freq_var_.set_num(num_freqs_);
    } else {
      power_var_.set_num(0);
      real_var_.set_num(0);
      imag_var_.set_num(0);
      freq_var_.set_num(0);
    }

#ifdef DEBUG
    cerr << "PowerResult::init(), computing power through a surface which is "
         << x_size_ << " x " << y_size_ << " x " << z_size_ 
         << " in size." << endl;
    cerr << "Frequency range: " << freq_start_ << " to " 
         << freq_stop_ << ", spacing: " << freq_space_ << ", number: "
         << num_freqs_ << endl;
    cerr << "Cell area is " << cell_area_ << endl;
#endif 
  } else {
#ifdef DEBUG
    cerr << "PowerResult::init(), rank " << MPI_RANK << " has no data to " 
         << "contribute. " << endl;
#endif
    power_var_.set_num(0);
    real_var_.set_num(0);
    imag_var_.set_num(0);
    freq_var_.set_num(0);
  }
}

void PowerResult::deinit()
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



  // Rearrange this so that there is only one loop traversing the
  // data, and put the time and freq calculations inside?

  // Compute the instantaneous power through the surface at this
  // instant in time.
  unsigned int idx = 0;
  time_power_ = 0;

  if (has_data_)
  {
    for (int x = xmin_; x < xmax_; x++) 
    {
      for (int y = ymin_; y < ymax_; y++) 
      {
        for (int z = zmin_; z < zmax_; z++) 
        {
          et1 = plane_->get_e_t1(x, y, z);
          et2 = plane_->get_e_t2(x, y, z);
        
          ht1 = (plane_->get_h_t1(x, y, z) 
                 + plane_->get_h_t1(x + step_x_, y + step_y_, 
                                    z + step_z_)) / 2;
          ht2 = (plane_->get_h_t2(x, y, z) 
                 + plane_->get_h_t2(x + step_x_, y + step_y_, 
                                    z + step_z_)) / 2;

          time_power_ += (et1 * ht2 - et2 * ht1) * cell_area_;
        }
      }
    }
  }

  MPI_Reduce(&time_power_, &time_power_, 1, GRID_MPI_TYPE, MPI_SUM, 0, 
             MPI_COMM_WORLD);

  // Compute the power throught the surface in the frequency domain
  for (unsigned int i = 0; i <= num_freqs_; i++)
  {
    cos_temp = cos(2 * PI * freqs_[i] * time);
    sin_temp = sin(2 * PI * freqs_[i] * time);

    p_real2 = 0;
    p_imag2 = 0;

    power_real_[i] = time_power_ * cos_temp;
    power_imag_[i] = time_power_ * sin_temp;
    

    // What was I thinking here???!?!?

//     unsigned int idx = 0;

//     if (has_data_)
//     {
//       for (int x = xmin_; x < xmax_; x++) 
//       {
//         for (int y = ymin_; y < ymax_; y++) 
//         {
//           for (int z = zmin_; z < zmax_; z++) 
//           {
//             et1 = plane_->get_e_t1(x, y, z);
//             et2 = plane_->get_e_t2(x, y, z);

//             ht1 = (plane_->get_h_t1(x, y, z) 
//                    + plane_->get_h_t1(x + step_x_, y + step_y_, 
//                                       z + step_z_)) / 2;
//             ht2 = (plane_->get_h_t2(x, y, z) 
//                    + plane_->get_h_t2(x + step_x_, y + step_y_, 
//                                       z + step_z_)) / 2;
          
//             et1_real = et1 * cos_temp;
//             et1_imag = (-1) * et1 * sin_temp;
          
//             et2_real = et2 * cos_temp;
//             et2_imag = (-1) * et2 * sin_temp;
          
//             ht1_real = ht1 * cos_temp;
//             ht1_imag = (-1) * ht1 * sin_temp;
          
//             ht2_real = ht2 * cos_temp;
//             ht2_imag = (-1) * ht2 * sin_temp;
          
//             et1r_[idx] += et1_real;
//             et2r_[idx] += et2_real;
//             et1i_[idx] += et1_imag;
//             et2i_[idx] += et2_imag;

//             ht1r_[idx] += ht1_real;
//             ht2r_[idx] += ht2_real;
//             ht1i_[idx] += ht1_imag;
//             ht2i_[idx] += ht2_imag;

//             p_real2 += (et1r_[idx] * ht2r_[idx] + et1i_[idx] * ht2i_[idx])
//               - (et2r_[idx] * ht1r_[idx] + et2i_[idx] * ht1i_[idx]);

//             p_imag2 += et1i_[idx] * ht2r_[idx] - ht2i_[idx] * et1r_[idx] 
//               + ht1i_[idx] * et2r_[idx] - et2i_[idx] * ht1r_[idx];

//             ++idx;
//           }
//         }
//       }
//     } 

//     MPI_Reduce(&p_real2, &p_real2, 1, GRID_MPI_TYPE, MPI_SUM, 0, 
//                MPI_COMM_WORLD);

//     MPI_Reduce(&p_imag2, &p_imag2, 1, GRID_MPI_TYPE, MPI_SUM, 0, 
//                MPI_COMM_WORLD);
    
//     power_real_[i] = 0.5 * p_real2 * cell_area_;
//     power_imag_[i] = 0.5 * p_imag2 * cell_area_;

  }

  return variables_;
}
