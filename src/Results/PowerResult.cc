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
#include "../PlaneTiling.hh"

#include <string.h> // for memset
#include <math.h>

PowerResult::PowerResult()
  : power_real_(0), power_imag_(0), time_power_(0),
    has_data_(false), xmin_(0), ymin_(0), zmin_(0), 
    xmax_(0), ymax_(0), zmax_(0),
    step_x_(0), step_y_(0), step_z_(0), 
    plane_(0), normal_(X_AXIS), 
    export_dfts_(false)
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
    power_real_(0), power_imag_(0), time_power_(0),
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

  /* For setting up dimensions when exporting DFT's */ 
  unsigned int global_t1_len, global_t2_len;
  unsigned int local_t1_len, local_t2_len;
  unsigned int t1_start, t2_start;

  /* Region must be in out local sub-domain */ 
  region_ = grid.get_local_region(*box_);
  shared_ptr<Block> global_b_ = grid.get_global_region(*box_);

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

    local_t1_len = y_size_;
    global_t1_len = (*global_b_).ylen();
    t1_start = (*region_).ystart() - (*global_b_).ystart();

    local_t2_len = z_size_;
    global_t2_len = (*global_b_).zlen();
    t2_start = (*region_).zstart() - (*global_b_).zstart();
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

    local_t1_len = z_size_;
    global_t1_len = (*global_b_).zlen();
    t1_start = (*region_).zstart() - (*global_b_).zstart();

    local_t2_len = x_size_;
    global_t2_len = (*global_b_).xlen();
    t2_start = (*region_).xstart() - (*global_b_).xstart();
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

    local_t1_len = x_size_;
    global_t1_len = (*global_b_).xlen();
    t1_start = (*region_).xstart() - (*global_b_).xstart();

    local_t2_len = y_size_;
    global_t2_len = (*global_b_).ylen();
    t2_start = (*region_).ystart() - (*global_b_).ystart();

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

    et1r_ = new field_t[sz];
    et1i_ = new field_t[sz];
    ht1r_ = new field_t[sz];
    ht1i_ = new field_t[sz];

    et2r_ = new field_t[sz];
    et2i_ = new field_t[sz];
    ht2r_ = new field_t[sz];
    ht2i_ = new field_t[sz];

    // new throws its own exception
//     if (!et1r_ || !et1i_ || !ht1r_ || !ht1i_
//         || !et2r_ || !et2i_ || !ht2r_ || !ht2i_)
//       throw MemoryException();

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
    power_real_ = new field_t[frequencies_.length()];
    power_imag_ = new field_t[frequencies_.length()];

//     if (!freqs_ || !power_real_ || !power_imag_)
//       throw MemoryException();

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

#ifdef DEBUG
    cerr << "PowerResult::init(), computing power through a surface which is "
         << x_size_ << " x " << y_size_ << " x " << z_size_ 
         << " in size." << endl;
    cerr << "Frequency range: " << frequencies_.get_start() << " to " 
         << frequencies_.get_end() << ", spacing: " 
         << frequencies_.get_spacing() << ", number: "
         << frequencies_.length() << endl;
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


  // Set up the optional DFT export
  if (export_dfts_)
  {
    variables_["et1_dft_r"] = &dfts_[0];
    variables_["et1_dft_i"] = &dfts_[1];
    variables_["et2_dft_r"] = &dfts_[2];
    variables_["et2_dft_i"] = &dfts_[3];
    variables_["ht1_dft_r"] = &dfts_[4];
    variables_["ht1_dft_i"] = &dfts_[5];
    variables_["ht2_dft_r"] = &dfts_[6];
    variables_["ht2_dft_i"] = &dfts_[7];

    for (int i = 0; i < 8; i++)
    {
      dfts_[i].reset();

      if (has_data_)
        dfts_[i].set_num(frequencies_.length() * x_size_ * y_size_ * z_size_);
      else
        dfts_[i].set_num(0);

      dfts_[i].has_time_dimension(false);
      dfts_[i].add_dimension("f", frequencies_.length(), frequencies_.length(), 0);
      dfts_[i].add_dimension("x", local_t1_len, global_t2_len, t1_start);
      dfts_[i].add_dimension("y", local_t2_len, global_t2_len, t2_start);
    }

    dfts_[0].set_name(base_name_ + "et1_dft_r");
    dfts_[1].set_name(base_name_ + "et1_dft_i");
    dfts_[2].set_name(base_name_ + "et2_dft_r");
    dfts_[3].set_name(base_name_ + "et2_dft_i");
    dfts_[4].set_name(base_name_ + "ht1_dft_r");
    dfts_[5].set_name(base_name_ + "ht1_dft_i");
    dfts_[6].set_name(base_name_ + "ht2_dft_r");
    dfts_[7].set_name(base_name_ + "ht2_dft_i");

    dfts_[0].set_ptr(et1r_);
    dfts_[1].set_ptr(et1i_);
    dfts_[2].set_ptr(et2r_);
    dfts_[3].set_ptr(et2i_);
    dfts_[4].set_ptr(ht1r_);
    dfts_[5].set_ptr(ht1i_);
    dfts_[6].set_ptr(ht2r_);
    dfts_[7].set_ptr(ht2i_);

  }

}

void PowerResult::deinit()
{
  if (plane_)
  {
    delete plane_;
    plane_ = 0;
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
  delta_t e_time = dt * time_step;
  delta_t h_time = dt * (static_cast<delta_t>(time_step) - 0.5);

  field_t e_cos_temp, e_sin_temp;
  field_t h_cos_temp, h_sin_temp;

  // Compute the instantaneous power through the surface at this
  // instant in time.
  unsigned int idx = 0;
  time_power_ = 0;

  if (has_data_)
  {
    // Make this a templated function taking the GridPlane as a template
    // parameter for speed, so we can avoid virtual function calls. 
    for (int x = xmin_; x < xmax_; x++) 
    {
      for (int y = ymin_; y < ymax_; y++) 
      {
        for (int z = zmin_; z < zmax_; z++) 
        {
          et1 = plane_->get_avg_e_t1(x, y, z);
          et2 = plane_->get_avg_e_t2(x, y, z);
        
          ht1 = plane_->get_avg_h_t1(x, y, z);
          ht2 = plane_->get_avg_h_t2(x, y, z);

          time_power_ += (et1 * ht2 - et2 * ht1) * cell_area_;
        }
      }
    }
  }

  MPI_Reduce(&time_power_, &time_power_, 1, GRID_MPI_TYPE, MPI_SUM, 0, 
             MPI_COMM_WORLD);

  // Compute the power throught the surface in the frequency domain
  for (unsigned int i = 0; i < frequencies_.length(); i++)
  {
    e_cos_temp = cos(-2 * PI * frequencies_.get(i) * e_time);
    e_sin_temp = sin(-2 * PI * frequencies_.get(i) * e_time);

    h_cos_temp = cos(-2 * PI * frequencies_.get(i) * h_time);
    h_sin_temp = sin(-2 * PI * frequencies_.get(i) * h_time);

    p_real2 = 0;
    p_imag2 = 0;

    unsigned int idx = i * x_size_ * y_size_ * z_size_;

    if (has_data_)
    {
      // Template this
      for (int x = xmin_; x < xmax_; x++) 
      {
        for (int y = ymin_; y < ymax_; y++) 
        {
          for (int z = zmin_; z < zmax_; z++) 
          {
            et1 = plane_->get_avg_e_t1(x, y, z);
            et2 = plane_->get_avg_e_t2(x, y, z);
            
            ht1 = plane_->get_avg_h_t1(x, y, z);
            ht2 = plane_->get_avg_h_t2(x, y, z);

            // Replace with complex numbers?
            et1_real = et1 * e_cos_temp;
            et1_imag = (-1) * et1 * e_sin_temp;
          
            et2_real = et2 * e_cos_temp;
            et2_imag = (-1) * et2 * e_sin_temp;

            ht1_real = ht1 * h_cos_temp;
            ht1_imag = (-1) * ht1 * h_sin_temp;

            ht2_real = ht2 * h_cos_temp;
            ht2_imag = (-1) * ht2 * h_sin_temp;
          
            et1r_[idx] += et1_real;
            et2r_[idx] += et2_real;
            et1i_[idx] += et1_imag;
            et2i_[idx] += et2_imag;

            ht1r_[idx] += ht1_real;
            ht2r_[idx] += ht2_real;
            ht1i_[idx] += ht1_imag;
            ht2i_[idx] += ht2_imag;

            // Hmmm....
            p_real2 += ((et1r_[idx] * ht2r_[idx] + et1i_[idx] * ht2i_[idx])
              - (et2r_[idx] * ht1r_[idx] + et2i_[idx] * ht1i_[idx])) 
              * cell_area_;

            p_imag2 += ((et1i_[idx] * ht2r_[idx] - ht2i_[idx] * et1r_[idx]) 
              + (ht1i_[idx] * et2r_[idx] - et2i_[idx] * ht1r_[idx]))
              * cell_area_;

            ++idx;
          }
        }
      }
    } 

    MPI_Reduce(&p_real2, &p_real2, 1, GRID_MPI_TYPE, MPI_SUM, 0, 
               MPI_COMM_WORLD);

    MPI_Reduce(&p_imag2, &p_imag2, 1, GRID_MPI_TYPE, MPI_SUM, 0, 
               MPI_COMM_WORLD);
    
    power_real_[i] = p_real2;
    power_imag_[i] = p_imag2;

  }

  return variables_;
}

ostream& PowerResult::to_string(ostream &os) const
{
  os << "PowerResult returning power passing through the "
     << face_string(face_) << " of a "; // << box_->get();

  return os;
}

class TimePowerAlg
{
public:
  class Data
  {
  public:
    field_t area_;
    field_t tp_; 

    Data(field_t area)
      : area_(area), tp_(0)
    {}
  }

  static inline alg(Fields_t &f, Data &data)
  {
    data.tp_ += (f.et1_avg * f.ht2_avg - f.et2_avg * f.ht1_avg) 
      * data.area_;
  }
};

class DFTPowerAlg
{
public:
  class Data
  {
  public:
    int idx;
    
    field_t e_cos_temp;
    field_t e_sin_temp;

    field_t h_cos_temp;
    field_t h_sin_temp;

    field_t *et1r_, *et2r_, *et1i_, *et2i_;
    field_t *ht1r_, *ht2r_, *ht1i_, *ht2i_;

    field_t p_real;
    field_t p_imag;

    field_t cell_area;
  };

  static inline void alg(Fields_t &f, Data &data)
  {
    data.et1r_[data.idx] += f.et1_avg * data.e_cos_temp;
    data.et1i_[data.idx] += -1 * f.et1_avg * data.e_sin_temp;

    data.et2r_[data.idx] += f.et2_avg * data.e_cos_temp;
    data.et2i_[data.idx] += -1 * f.et2_avg * data.e_sin_temp;

    data.ht1r_[data.idx] += f.ht1_avg * data.h_cos_temp;
    data.ht1i_[data.idx] += -1 * f.ht1_avg * data.h_sin_temp;

    data.ht2r_[data.idx] += f.ht2_avg * data.h_cos_temp;
    data.ht2i_[data.idx] += -1 * f.ht2_avg * data.h_sin_temp;

    data.p_real += ((data.et1r_[data.idx] * data.ht2r_[data.idx] 
                     + data.et1i_[data.idx] * data.ht2i_[data.idx])
                - (data.et2r_[data.idx] * data.ht1r_[data.idx] 
                   + data.et2i_[data.idx] * data.ht1i_[data.idx])) 
      * data.cell_area;

    data.p_imag += ((data.et1i_[data.idx] * data.ht2r_[data.idx] 
                     - data.ht2i_[data.idx] * data.et1r_[data.idx]) 
                + (data.ht1i_[data.idx] * data.et2r_[data.idx] 
                   - data.et2i_[data.idx] * data.ht1r_[data.idx]))
      * data.cell_area;

    ++data.idx;
  }
};

// An experimental implementation of get_result which uses algorithm
// objects and PlaneTiling.
map<string, Variable *> &PowerResult::get_result_test(const Grid &grid, 
                                                      unsigned int time_step)
{
  time_power_ = 0;
  field_t p_real = 0;
  field_t p_imag = 0;

  if (has_data_)
  {
    // Make this a templated function taking the GridPlane as a template
    // parameter for speed, so we can avoid virtual function calls. 

    TimePowerAlg::Data data(cell_area_);
    
    PlaneTiling<TimePowerAlg, TimePowerAlg::Data>::loop(grid, (*region_),
                                                        face_, data);
    time_power_ = data.tp_;
  }

  MPI_Reduce(&time_power_, &time_power_, 1, GRID_MPI_TYPE, MPI_SUM, 0, 
             MPI_COMM_WORLD);

  DFTPowerAlg::Data dftdata;

  dftdata.cell_area = cell_area_;

  delta_t dt = grid.get_deltat();
  delta_t e_time = dt * time_step;
  delta_t h_time = dt * (static_cast<delta_t>(time_step) - 0.5);

  for (unsigned int i = 0; i < frequencies_.length(); i++)
  {
    dftdata.p_real = 0;
    dftdata.p_imag = 0;

    if (has_data_)
    {
      dftdata.e_cos_temp = cos(-2 * PI * frequencies_.get(i) * e_time);
      dftdata.e_sin_temp = sin(-2 * PI * frequencies_.get(i) * e_time);

      dftdata.h_cos_temp = cos(-2 * PI * frequencies_.get(i) * h_time);
      dftdata.h_sin_temp = sin(-2 * PI * frequencies_.get(i) * h_time);

      dftdata.p_real = 0;
      dftdata.p_imag = 0;

      dftdata.idx = i * x_size_ * y_size_ * z_size_;
      
      PlaneTiling<DFTPowerAlg, DFTPowerAlg::Data>::loop(grid, (*region_),
                                                        face_, dftdata);
    }
    
    MPI_Reduce(&dftdata.p_real, &dftdata.p_real, 1, 
               GRID_MPI_TYPE, MPI_SUM, 0, 
               MPI_COMM_WORLD);
    
    MPI_Reduce(&dftdata.p_imag, &dftdata.p_imag, 1, 
               GRID_MPI_TYPE, MPI_SUM, 0, 
               MPI_COMM_WORLD);

    power_real_[i] = dftdata.p_real;
    power_imag_[i] = dftdata.p_imag;
  }

  return variables_;

}
