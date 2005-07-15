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

#include "AvgPowerResult.hh"
#include "../Globals.hh"
#include "../Constants.hh"
#include "../PlaneTiling.hh"

#include <cstring>
#include <cmath>

AvgPowerResult::AvgPowerResult()
  : power_real_(0), power_imag_(0), 
    et1_(0), et2_(0), ht1_(0), ht2_(0),
    prev_et1_(0), prev_et2_(0), 
    cos_temp_(0), sin_temp_(0),
    e_cos_temp_(0), e_sin_temp_(0)
{}

AvgPowerResult::AvgPowerResult(field_t freq_start, field_t freq_stop, 
                               unsigned int num_freqs)
  : DFTResult(freq_start, freq_stop, num_freqs), 
    power_real_(0), power_imag_(0), 
    et1_(0), et2_(0), ht1_(0), ht2_(0),
    prev_et1_(0), prev_et2_(0),
    cos_temp_(0), sin_temp_(0),
    e_cos_temp_(0), e_sin_temp_(0)
{}

AvgPowerResult::~AvgPowerResult()
{
  deinit();
}

void AvgPowerResult::init(const Grid &grid)
{
  real_var_.reset();
  imag_var_.reset();
  freq_var_.reset();

  real_var_.has_time_dimension(false);
  imag_var_.has_time_dimension(false);
  freq_var_.has_time_dimension(false);

  post_vars_["real_power"] = &real_var_;
  post_vars_["imag_power"] = &imag_var_;
  pre_vars_["freqs"] = &freq_var_;

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

  real_var_.set_name(base_name_ + "_real");
  imag_var_.set_name(base_name_ + "_imag");
  freq_var_.set_name(base_name_ + "_freqs");
  
  real_var_.add_dimension("Frequency", frequencies_.length(), 
                          frequencies_.length(), 0);
  imag_var_.add_dimension("Frequency", frequencies_.length(), 
                          frequencies_.length(), 0);
  freq_var_.add_dimension("Frequency", frequencies_.length(), 
                          frequencies_.length(), 0);
  
  if (!setup_only && has_data_)
  {
    cos_temp_ = new field_t[frequencies_.length()];
    sin_temp_ = new field_t[frequencies_.length()];

    e_cos_temp_ = new field_t[frequencies_.length()];
    e_sin_temp_ = new field_t[frequencies_.length()];

    memset(cos_temp_, 0, sizeof(field_t) * frequencies_.length());
    memset(sin_temp_, 0, sizeof(field_t) * frequencies_.length());

    memset(e_cos_temp_, 0, sizeof(field_t) * frequencies_.length());
    memset(e_sin_temp_, 0, sizeof(field_t) * frequencies_.length());

    /* GRR, ARGH! */
    unsigned int sz = frequencies_.length() * x_size_ * y_size_ * z_size_;

    et1_ = new complex<field_t>[sz];
    et2_ = new complex<field_t>[sz];

    ht1_ = new complex<field_t>[sz];
    ht2_ = new complex<field_t>[sz];

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

    if (MPI_RANK == 0) // All data is collected to rank 0 by this
                       // result, rather than by the DataWriters.
    {
      real_var_.set_num(frequencies_.length());
      imag_var_.set_num(frequencies_.length());
      freq_var_.set_num(frequencies_.length());
    } else {
      real_var_.set_num(0);
      imag_var_.set_num(0);
      freq_var_.set_num(0);
    }

  } else {
    real_var_.set_num(0);
    imag_var_.set_num(0);
    freq_var_.set_num(0);
  }
}

void AvgPowerResult::deinit()
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

  if (et1_)
  {
    delete[] et1_;
    delete[] et2_;
    delete[] ht1_;
    delete[] ht2_;

    et1_ = 0;
    et2_ = 0;
    ht1_ = 0;
    ht2_ = 0;
  }

  if (prev_et1_)
  {
    delete[] prev_et1_;
    prev_et1_ = 0;

    delete[] prev_et2_;
    prev_et2_ = 0;
  }

  if (cos_temp_)
  {
    delete[] cos_temp_;
    delete[] sin_temp_;

    delete[] e_cos_temp_;
    delete[] e_sin_temp_;

    cos_temp_ = 0;
    sin_temp_ = 0;
    e_cos_temp_ = 0;
    e_sin_temp_ = 0;
  }
}

/**
 * This class is a meta-program template which calculates the running
 * DFT of the E and H fields of interest.
 */ 
class DFTPowerAlg
{
public:
  class Data
  {
  public:
    int idx;
    
    field_t cos_temp;
    field_t sin_temp;

    field_t e_cos_temp;
    field_t e_sin_temp;

    complex<field_t> *et1_, *et2_, *ht1_, *ht2_;

    const field_t *prev_et1_, *prev_et2_;

  };

  static inline void alg(const int &x, const int &y, const int &z, 
                         Fields_t &f, Data &data)
  {
    // It is necessary to store the old value of E so that the average
    // of two time steps can be calculated, otherwise the result will
    // be incorrect.
    field_t et1_tavg = (f.et1_avg); // + data.prev_et1_[data.idx]) / 2;
    field_t et2_tavg = (f.et2_avg); // + data.prev_et2_[data.idx]) / 2;

    // THE TIME AVERAGING SEEMS TO BE CAUSING PROBLEMS
    data.et1_[data.idx] += complex<field_t>(et1_tavg * data.e_cos_temp, 
                                            -1 * et1_tavg * data.e_sin_temp);
    data.et2_[data.idx] += complex<field_t>(et2_tavg * data.e_cos_temp, 
                                            -1 * et2_tavg * data.e_sin_temp);
    data.ht1_[data.idx] += complex<field_t>(f.ht1_avg * data.cos_temp,
                                            -1 * f.ht1_avg * data.sin_temp);
    data.ht2_[data.idx] += complex<field_t>(f.ht2_avg * data.cos_temp, 
                                            -1 * f.ht2_avg * data.sin_temp);

    ++data.idx;
  }
};

/**
 * This class is a meta-program template which calculates the running
 * DFT of the E and H fields of interest.
 *
 * This version is hopefully faster...
 */ 
class DFTPowerAlgFast
{
public:
  class Data
  {
  public:
    int idx;
    
    field_t *cos_temp;
    field_t *sin_temp;

    field_t *e_cos_temp;
    field_t *e_sin_temp;

    complex<field_t> *et1_, *et2_, *ht1_, *ht2_;

    const field_t *prev_et1_, *prev_et2_;

    int num_f_;
  };

  static inline void alg(const int &x, const int &y, const int &z, 
                         Fields_t &f, Data &data)
  {
    // It is necessary to store the old value of E so that the average
    // of two time steps can be calculated, otherwise the result will
    // be incorrect.
    field_t et1_tavg = (f.et1_avg); // + data.prev_et1_[data.idx]) / 2;
    field_t et2_tavg = (f.et2_avg); // + data.prev_et2_[data.idx]) / 2;

    for (unsigned int i = 0; i < data.num_f_; i++)
    {
      // THE TIME AVERAGING SEEMS TO BE CAUSING PROBLEMS
      data.et1_[data.idx] += complex<field_t>
        (et1_tavg * data.e_cos_temp[i], 
         -1 * et1_tavg * data.e_sin_temp[i]);
      data.et2_[data.idx] += complex<field_t>
        (et2_tavg * data.e_cos_temp[i], 
         -1 * et2_tavg * data.e_sin_temp[i]);
      data.ht1_[data.idx] += complex<field_t>
        (f.ht1_avg * data.cos_temp[i],
         -1 * f.ht1_avg * data.sin_temp[i]);
      data.ht2_[data.idx] += complex<field_t>
        (f.ht2_avg * data.cos_temp[i], 
         -1 * f.ht2_avg * data.sin_temp[i]);

      ++data.idx;
    }
  }
};

/**
 * This class calculates the total power through the plane after all
 * the time steps have been completed.
 */
struct CalcPowerAlg
{
  struct Data
  {
    int idx;
    int stride;

    complex<field_t> *et1_, *et2_, *ht1_, *ht2_;

    field_t p_real;
    field_t p_imag;

    field_t cell_area;
  };

  static inline void alg(const int &x, const int &y, const int &z,
                         Fields_t &f, Data &data)
  {
    complex<field_t> temp = (data.et1_[data.idx] * conj(data.ht2_[data.idx])
      - data.et2_[data.idx] * conj(data.ht1_[data.idx])) * data.cell_area;
 
    // Cheap, approximate integration.
    data.p_real += temp.real();
    data.p_imag += temp.imag();


    data.idx += data.stride;
  }
};

/**
 * This class updates the stored value of the E tangential components 
 */
// class PrevEupdate
// {
// public:
//   class Data 
//   {
//   public:
//     field_t *prev_et1_;
//     field_t *prev_et2_;
    
//     int idx;

//     Data(field_t *pet1, field_t *pet2)
//       : prev_et1_(pet1), prev_et2_(pet2), idx(0)
//     {}
//   };

//   static inline void alg(const int &x, const int &y, const int &z,
//                          Fields_t &f, Data &data)
//   {
//     data.prev_et1_[data.idx] = f.et1_avg;
//     data.prev_et2_[data.idx] = f.et2_avg;

//     ++data.idx;
//   }
// };

void AvgPowerResult::calculate_result(const Grid &grid, 
                                      unsigned int time_step)
{
//   DFTPowerAlg::Data dftdata;

//   dftdata.et1_ = et1_;
//   dftdata.et2_ = et2_;

//   dftdata.ht1_ = ht1_;
//   dftdata.ht2_ = ht2_;

//   dftdata.prev_et1_ = prev_et1_;
//   dftdata.prev_et2_ = prev_et2_;

//   delta_t dt = grid.get_deltat();
//   delta_t time = dt * time_step;
//   delta_t e_time = dt * (static_cast<delta_t>(time_step) + 0.5);

//   for (unsigned int i = 0; i < frequencies_.length(); i++)
//   {
//     if (has_data_)
//     {
//       dftdata.cos_temp = cos(-2 * PI * frequencies_.get(i) * time);
//       dftdata.sin_temp = sin(-2 * PI * frequencies_.get(i) * time);

//       dftdata.e_cos_temp = cos(-2 * PI * frequencies_.get(i) * e_time);
//       dftdata.e_sin_temp = sin(-2 * PI * frequencies_.get(i) * e_time);

//       dftdata.idx = i * x_size_ * y_size_ * z_size_;
      
//       PlaneTiling<DFTPowerAlg, DFTPowerAlg::Data>::loop(grid, (*region_),
//                                                         face_, dftdata);
//     }
//   }

  DFTPowerAlgFast::Data dftdata;

  dftdata.et1_ = et1_;
  dftdata.et2_ = et2_;

  dftdata.ht1_ = ht1_;
  dftdata.ht2_ = ht2_;

  dftdata.prev_et1_ = prev_et1_;
  dftdata.prev_et2_ = prev_et2_;

  dftdata.cos_temp = cos_temp_;
  dftdata.sin_temp = sin_temp_;

  dftdata.e_cos_temp = e_cos_temp_;
  dftdata.e_sin_temp = e_sin_temp_;

  dftdata.idx = 0;

  dftdata.num_f_ = frequencies_.length();

  delta_t dt = grid.get_deltat();
  delta_t time = dt * time_step;
  delta_t e_time = dt * (static_cast<delta_t>(time_step) + 0.5);

  for (unsigned int i = 0; i < frequencies_.length(); i++)
  {
    if (has_data_)
    {
      dftdata.cos_temp[i] = cos(-2 * PI * frequencies_.get(i) * time);
      dftdata.sin_temp[i] = sin(-2 * PI * frequencies_.get(i) * time);

      dftdata.e_cos_temp[i] = cos(-2 * PI * frequencies_.get(i) * e_time);
      dftdata.e_sin_temp[i] = sin(-2 * PI * frequencies_.get(i) * e_time);
    }
  }

  if (has_data_)
    PlaneTiling<DFTPowerAlgFast, DFTPowerAlgFast::Data>
      ::loop(grid, (*region_), face_, dftdata);
  
  // Store the current value of Et1, Et2
//   PrevEupdate::Data pedata(prev_et1_, prev_et2_);
//   PlaneTiling<PrevEupdate, PrevEupdate::Data>::loop(grid, (*region_),
//                                                     face_, pedata);
}

void AvgPowerResult::calculate_post_result(const Grid &grid)
{
  CalcPowerAlg::Data data;

  data.et1_ = et1_;
  data.et2_ = et2_;

  data.ht1_ = ht1_;
  data.ht2_ = ht2_;

  data.cell_area = cell_area_;

  for (unsigned int i = 0; i < frequencies_.length(); i++)
  {
    data.p_real = 0;
    data.p_imag = 0;

    if (has_data_)
    {
      data.idx = i; 
      data.stride = frequencies_.length(); 

      PlaneTiling<CalcPowerAlg, CalcPowerAlg::Data>
        ::loop(grid, (*region_), face_, data);
    }

    field_t temp_real, temp_imag;
    MPI_Reduce(&data.p_real, &temp_real, 1, 
               GRID_MPI_TYPE, MPI_SUM, 0, 
               MPI_COMM_PHRED);
    data.p_real = temp_real;

    MPI_Reduce(&data.p_imag, &temp_imag, 1, 
               GRID_MPI_TYPE, MPI_SUM, 0, 
               MPI_COMM_PHRED);
    data.p_imag = temp_imag;

    power_real_[i] = data.p_real;
    power_imag_[i] = data.p_imag;
  }
}

ostream& AvgPowerResult::to_string(ostream &os) const
{
  os << "AvgPowerResult returning power passing through the "
     << face_string(face_) << " of a box in grid coordinates " 
     << *region_ << endl;

  return os;
}

