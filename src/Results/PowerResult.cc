/* 
   Phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004-2005 Matt Hughes <mhughe@uvic.ca>

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
#include "../PlaneTiling.hh"

#include <cstring>
#include <cmath>

PowerResult::PowerResult()
  : power_real_(0), power_imag_(0), time_power_(0),
  has_data_(false) //, et1_(0), et2_(0), ht1_(0), ht2_(0)
{
  real_var_.has_time_dimension(false);
  imag_var_.has_time_dimension(false);
  freq_var_.has_time_dimension(false);
  power_var_.has_time_dimension(true);

  post_vars_["real_power"] = &real_var_;
  post_vars_["imag_power"] = &imag_var_;
  pre_vars_["freqs"] = &freq_var_;
  variables_["time_power"] = &power_var_;
}

PowerResult::PowerResult(field_t freq_start, field_t freq_stop, 
                         unsigned int num_freqs)
  : DFTResult(freq_start, freq_stop, num_freqs), 
    power_real_(0), power_imag_(0), time_power_(0),
  has_data_(false) //, et1_(0), et2_(0), ht1_(0), ht2_(0)
{
  real_var_.has_time_dimension(false);
  imag_var_.has_time_dimension(false);
  freq_var_.has_time_dimension(false);
  power_var_.has_time_dimension(true);

  post_vars_["real_power"] = &real_var_;
  post_vars_["imag_power"] = &imag_var_;
  pre_vars_["freqs"] = &freq_var_;
  variables_["time_power"] = &power_var_;
}

PowerResult::~PowerResult()
{
  deinit();
}

void PowerResult::init(const Grid &grid)
{
  /* Region must be in out local sub-domain */ 
  shared_ptr<CellSet> cells = grid.get_cellset(*box_);

  region_ = cells->get_local_block();
  shared_ptr<Block> global_b = cells->get_global_block();

  x_size_ = (*region_).xlen();
  y_size_ = (*region_).ylen();
  z_size_ = (*region_).zlen();

  /* Region must be a plane; set up grid plane */
  if (face_ == FRONT || face_ == BACK)
  {
    x_size_ = 1;
    cell_area_ = grid.get_deltay() * grid.get_deltaz();

    if (face_ == FRONT && (*region_).has_face_data(FRONT))
    {
      has_data_ = true;
    }

    if (face_ == BACK && (*region_).has_face_data(BACK))
    {
      has_data_ = true;
    }
  } 
  else if (face_ == LEFT || face_ == RIGHT) 
  {
    y_size_ = 1;
    cell_area_ = grid.get_deltax() * grid.get_deltaz();

    if (face_ == RIGHT && (*region_).has_face_data(RIGHT))
    {
      has_data_ = true;
    }

    if (face_ == LEFT && (*region_).has_face_data(LEFT))
    {
      has_data_ = true;
    }
  } 
  else if (face_ == TOP || face_ == BOTTOM)
  {
    z_size_ = 1;
    cell_area_ = grid.get_deltax() * grid.get_deltay();

    if (face_ == TOP && (*region_).has_face_data(TOP))
    {
      has_data_ = true;
    }

    if (face_ == BOTTOM && (*region_).has_face_data(BOTTOM))
    {
      has_data_ = true;
    }
  } else {
    throw ResultException("Region must be limited to a plane "
                          "for a PowerResult.");
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
  
  real_var_.add_dimension("Frequency", frequencies_.length(), 
                          frequencies_.length(), 0);
  imag_var_.add_dimension("Frequency", frequencies_.length(), 
                          frequencies_.length(), 0);
  freq_var_.add_dimension("Frequency", frequencies_.length(), 
                          frequencies_.length(), 0);
  power_var_.add_dimension("Power", 1, 1, 0);
  
  if (has_data_)
  {
    /* GRR, ARGH! */
    unsigned int sz = frequencies_.length() * x_size_ * y_size_ * z_size_;

//     et1_ = new complex<field_t>[sz];
//     et2_ = new complex<field_t>[sz];

//     ht1_ = new complex<field_t>[sz];
//     ht2_ = new complex<field_t>[sz];

    prev_et1_ = new field_t[sz];
    prev_et2_ = new field_t[sz];

    memset(prev_et1_, 0, sizeof(field_t) * sz);
    memset(prev_et2_, 0, sizeof(field_t) * sz);

    /* Set up the frequencies */ 
    power_real_ = new field_t[frequencies_.length()];
    power_imag_ = new field_t[frequencies_.length()];

    memset(power_imag_, 0, sizeof(field_t) * (frequencies_.length()));
    memset(power_real_, 0, sizeof(field_t) * (frequencies_.length()));

    /* Set up output variables. */
    imag_var_.set_ptr(power_imag_);
    real_var_.set_ptr(power_real_);  
    freq_var_.set_ptr(frequencies_.get_ptr());
    power_var_.set_ptr(&time_power_);

    if (MPI_RANK == 0) // All data is collected to rank 0 by this
                       // result, rather than by the DataWriters.
    {
      power_var_.set_num(1);
      real_var_.set_num(frequencies_.length());
      imag_var_.set_num(frequencies_.length());
      freq_var_.set_num(frequencies_.length());
    } else {
      power_var_.set_num(0);
      real_var_.set_num(0);
      imag_var_.set_num(0);
      freq_var_.set_num(0);
    }

// #ifdef DEBUG
//     cerr << "PowerResult::init(), computing power through a surface which is "
//          << x_size_ << " x " << y_size_ << " x " << z_size_ 
//          << " in size." << endl;
//     cerr << "Frequency range: " << frequencies_.get_start() << " to " 
//          << frequencies_.get_end() << ", spacing: " 
//          << frequencies_.get_spacing() << ", number: "
//          << frequencies_.length() << endl;
//     cerr << "Cell area is " << cell_area_ << endl;
// #endif 
  } else {
// #ifdef DEBUG
//     cerr << "PowerResult::init(), rank " << MPI_RANK << " has no data to " 
//          << "contribute. " << endl;
// #endif
    power_var_.set_num(0);
    real_var_.set_num(0);
    imag_var_.set_num(0);
    freq_var_.set_num(0);
  }
}

void PowerResult::deinit()
{
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

  if (prev_et1_)
  {
    delete[] prev_et1_;
    prev_et1_ = 0;

    delete[] prev_et2_;
    prev_et2_ = 0;
  }
}

ostream& PowerResult::to_string(ostream &os) const
{
  os << "PowerResult returning power passing through the "
     << face_string(face_) << " of a box in grid coordinates " 
     << *region_ << endl;

  return os;
}

/**
 * This class is a meta-program template which calculates the power
 * through plane in the time domain.
 */ 
class TimePowerAlg
{
public:
  class Data
  {
  public:
    field_t area_;
    field_t tp_; 

    const field_t *prev_et1_, *prev_et2_;

    int idx;

    Data(field_t area, const field_t *pet1, const field_t *pet2)
      : area_(area), tp_(0), prev_et1_(pet1), prev_et2_(pet2), idx(0)
    {}
  };

  static inline void alg(const int &x, const int &y, const int &z, 
                         Fields_t &f, Data &data)
  {
    // It is necessary to store the old value of E so that the average
    // of two time steps can be calculated, otherwise the result will
    // be incorrect.
    field_t et1_tavg = (f.et1_avg + data.prev_et1_[data.idx]) / 2;
    field_t et2_tavg = (f.et2_avg + data.prev_et2_[data.idx]) / 2;

    data.tp_ += (et1_tavg * f.ht2_avg - et2_tavg * f.ht1_avg) 
      * data.area_;

    ++data.idx;
  }
};

/**
 * This class updates the stored value of the E tangential components 
 */
class PrevEupdate
{
public:
  class Data 
  {
  public:
    field_t *prev_et1_;
    field_t *prev_et2_;
    
    int idx;

    Data(field_t *pet1, field_t *pet2)
      : prev_et1_(pet1), prev_et2_(pet2), idx(0)
    {}
  };

  static inline void alg(const int &x, const int &y, const int &z,
                         Fields_t &f, Data &data)
  {
    data.prev_et1_[data.idx] = f.et1_avg;
    data.prev_et2_[data.idx] = f.et2_avg;

    ++data.idx;
  }
};

void PowerResult::calculate_result(const Grid &grid, 
                                   unsigned int time_step)
{
  time_power_ = 0;

  if (has_data_)
  {
    TimePowerAlg::Data data(cell_area_, prev_et1_, prev_et2_);
    
    PlaneTiling<TimePowerAlg, TimePowerAlg::Data>::loop(grid, (*region_),
                                                        face_, data);
    time_power_ = data.tp_;
  }

  field_t temp_time_power;
  MPI_Reduce(&time_power_, &temp_time_power, 1, GRID_MPI_TYPE, MPI_SUM, 0, 
             MPI_COMM_PHRED);
  time_power_ = temp_time_power;

  // DFTPowerAlg::Data dftdata;

//   dftdata.cell_area = cell_area_;

  delta_t dt = grid.get_deltat();

  // I think these may be backwards
  delta_t e_time = dt * time_step;
  delta_t h_time = dt * (static_cast<delta_t>(time_step) - 0.5);

//   dftdata.et1_ = et1_;
//   dftdata.et2_ = et2_;

//   dftdata.ht1_ = ht1_;
//   dftdata.ht2_ = ht2_;

//   dftdata.prev_et1_ = prev_et1_;
//   dftdata.prev_et2_ = prev_et2_;

  for (unsigned int i = 0; i < frequencies_.length(); i++)
  {
//     dftdata.p_real = 0;
//     dftdata.p_imag = 0;

//     if (has_data_)
//     {
//       dftdata.e_cos_temp = cos(-2 * PI * frequencies_.get(i) * e_time);
//       dftdata.e_sin_temp = sin(-2 * PI * frequencies_.get(i) * e_time);

//       dftdata.h_cos_temp = cos(-2 * PI * frequencies_.get(i) * h_time);
//       dftdata.h_sin_temp = sin(-2 * PI * frequencies_.get(i) * h_time);

//       dftdata.p_real = 0;
//       dftdata.p_imag = 0;

//       dftdata.idx = i * x_size_ * y_size_ * z_size_;
      
//       PlaneTiling<DFTPowerAlg, DFTPowerAlg::Data>::loop(grid, (*region_),
//                                                         face_, dftdata);
//     }
    
//     field_t temp_real, temp_imag;
//     MPI_Reduce(&dftdata.p_real, &temp_real, 1, 
//                GRID_MPI_TYPE, MPI_SUM, 0, 
//                MPI_COMM_PHRED);
//     dftdata.p_real = temp_real;

//     MPI_Reduce(&dftdata.p_imag, &temp_imag, 1, 
//                GRID_MPI_TYPE, MPI_SUM, 0, 
//                MPI_COMM_PHRED);
//     dftdata.p_imag = temp_imag;

//     power_real_[i] = dftdata.p_real;
//     power_imag_[i] = dftdata.p_imag;

    power_real_[i] += time_power_ * cos(-2 * PI * frequencies_.get(i) * e_time);
    power_imag_[i] += time_power_ * sin(-2 * PI * frequencies_.get(i) * e_time);
  }

  // Store the current value of Et1, Et2
  PrevEupdate::Data pedata(prev_et1_, prev_et2_);
  PlaneTiling<PrevEupdate, PrevEupdate::Data>::loop(grid, (*region_),
                                                    face_, pedata);

}
