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

#include "SurfaceCurrentResult.hh"
#include "../GridPlane.hh"

#include <cmath>

SurfaceCurrentResult::SurfaceCurrentResult()
  : freq_start_(0), freq_stop_(0), num_freqs_(0), freq_space_(0), 
    freqs_(0), do_dft_(false)
{
  for (int i = 0; i < 6; i++)
  {
    Jt1_data_[i] = 0;
    Jt2_data_[i] = 0;
    Mt1_data_[i] = 0;
    Mt2_data_[i] = 0;
  }
}

SurfaceCurrentResult::~SurfaceCurrentResult()
{ deinit(); }

void SurfaceCurrentResult::init(const Grid &grid)
{
  variables_["back_Jy"] = &Jt1_[BACK];
  variables_["back_Jz"] = &Jt2_[BACK];
  variables_["back_My"] = &Mt1_[BACK];
  variables_["back_Mz"] = &Mt2_[BACK];

  variables_["front_Jy"] = &Jt1_[FRONT];
  variables_["front_Jz"] = &Jt2_[FRONT];
  variables_["front_My"] = &Mt1_[FRONT];
  variables_["front_Mz"] = &Mt2_[FRONT];

  variables_["left_Jx"] = &Jt1_[LEFT];
  variables_["left_Jz"] = &Jt2_[LEFT];
  variables_["left_Mx"] = &Mt1_[LEFT];
  variables_["left_Mz"] = &Mt2_[LEFT];

  variables_["right_Jx"] = &Jt1_[RIGHT];
  variables_["right_Jz"] = &Jt2_[RIGHT];
  variables_["right_Mx"] = &Mt1_[RIGHT];
  variables_["right_Mz"] = &Mt2_[RIGHT];

  variables_["bottom_Jx"] = &Jt1_[BOTTOM];
  variables_["bottom_Jy"] = &Jt2_[BOTTOM];
  variables_["bottom_Mx"] = &Mt1_[BOTTOM];
  variables_["bottom_My"] = &Mt2_[BOTTOM];

  variables_["top_Jx"] = &Jt1_[TOP];
  variables_["top_Jy"] = &Jt2_[TOP];
  variables_["top_Mx"] = &Mt1_[TOP];
  variables_["top_My"] = &Mt2_[TOP];

  if (box_.get())
  {
    region_ = grid.get_local_region(*(box_.get()));
  } else {
    throw ResultException("SurfaceCurrentResult has no surface defined!");
  }
  
  shared_ptr<Block> gregion = grid.get_global_region(*(box_.get()));

  if (do_dft_ && num_freqs_ > 2 && freq_stop_ > freq_start_)
  {
    freq_space_ = (freq_stop_ - freq_start_) / (num_freqs_ - 1);

    freqs_ = new field_t[num_freqs_];

    if (!freqs_)
      throw MemoryException();

    for (unsigned int i = 0; i < num_freqs_; i++)
    {
      freqs_[i] = freq_start_ + i * freq_space_;
    }
  }
  else 
    do_dft_ = false;
  
  unsigned int face_size = 0;

  for (int i = 0; i < 6; i++)
  {
    switch (i)
    {
    case FRONT:
      Jt1_[i].set_name(base_name_ + "_front_Jy");
      Jt2_[i].set_name(base_name_ + "_front_Jz");
      Mt1_[i].set_name(base_name_ + "_front_My");
      Mt2_[i].set_name(base_name_ + "_front_Mz");
      //face_size = (*region_).xmax() - (*region_).xmin();
      face_size = (*region_).ylen() * (*region_).zlen();

      Jt1_[i].add_dimension("y", (*region_).ylen(), (*gregion).ylen(),
                            (*region_).ystart());
      Jt1_[i].add_dimension("z", (*region_).zlen(), (*gregion).zlen(),
                            (*region_).zstart());

      Jt2_[i].add_dimension("y", (*region_).ylen(), (*gregion).ylen(),
                            (*region_).ystart());
      Jt2_[i].add_dimension("z", (*region_).zlen(), (*gregion).zlen(),
                            (*region_).zstart());

      Mt1_[i].add_dimension("y", (*region_).ylen(), (*gregion).ylen(),
                            (*region_).ystart());
      Mt1_[i].add_dimension("z", (*region_).zlen(), (*gregion).zlen(),
                            (*region_).zstart());

      Mt2_[i].add_dimension("y", (*region_).ylen(), (*gregion).ylen(),
                            (*region_).ystart());
      Mt2_[i].add_dimension("z", (*region_).zlen(), (*gregion).zlen(),
                            (*region_).zstart());
      break;

    case BACK:
      Jt1_[i].set_name(base_name_ + "_back_Jy");
      Jt2_[i].set_name(base_name_ + "_back_Jz");
      Mt1_[i].set_name(base_name_ + "_back_My");
      Mt2_[i].set_name(base_name_ + "_back_Mz");
      //face_size = (*region_).xmax() - (*region_).xmin();
      face_size = (*region_).ylen() * (*region_).zlen();

      Jt1_[i].add_dimension("y", (*region_).ylen(), (*gregion).ylen(),
                            (*region_).ystart());
      Jt1_[i].add_dimension("z", (*region_).zlen(), (*gregion).zlen(),
                            (*region_).zstart());

      Jt2_[i].add_dimension("y", (*region_).ylen(), (*gregion).ylen(),
                            (*region_).ystart());
      Jt2_[i].add_dimension("z", (*region_).zlen(), (*gregion).zlen(),
                            (*region_).zstart());

      Mt1_[i].add_dimension("y", (*region_).ylen(), (*gregion).ylen(),
                            (*region_).ystart());
      Mt1_[i].add_dimension("z", (*region_).zlen(), (*gregion).zlen(),
                            (*region_).zstart());

      Mt2_[i].add_dimension("y", (*region_).ylen(), (*gregion).ylen(),
                            (*region_).ystart());
      Mt2_[i].add_dimension("z", (*region_).zlen(), (*gregion).zlen(),
                            (*region_).zstart());
      break;

    case LEFT:
      Jt1_[i].set_name(base_name_ + "_left_Jx");
      Jt2_[i].set_name(base_name_ + "_left_Jz");
      Mt1_[i].set_name(base_name_ + "_left_Mx");
      Mt2_[i].set_name(base_name_ + "_left_Mz");
      //face_size = (*region_).ymax() - (*region_).ymin();
      face_size = (*region_).xlen() * (*region_).zlen();

      Jt1_[i].add_dimension("z", (*region_).zlen(), (*gregion).zlen(),
                            (*region_).zstart());
      Jt1_[i].add_dimension("x", (*region_).xlen(), (*gregion).xlen(),
                            (*region_).xstart());

      Jt2_[i].add_dimension("z", (*region_).zlen(), (*gregion).zlen(),
                            (*region_).zstart());
      Jt2_[i].add_dimension("x", (*region_).xlen(), (*gregion).xlen(),
                            (*region_).xstart());

      Mt1_[i].add_dimension("z", (*region_).zlen(), (*gregion).zlen(),
                            (*region_).zstart());
      Mt1_[i].add_dimension("x", (*region_).xlen(), (*gregion).xlen(),
                            (*region_).xstart());

      Mt2_[i].add_dimension("z", (*region_).zlen(), (*gregion).zlen(),
                            (*region_).zstart());
      Mt2_[i].add_dimension("x", (*region_).xlen(), (*gregion).xlen(),
                            (*region_).xstart());
      break;

    case RIGHT: 
      Jt1_[i].set_name(base_name_ + "_right_Jx");
      Jt2_[i].set_name(base_name_ + "_right_Jz");
      Mt1_[i].set_name(base_name_ + "_right_Mx");
      Mt2_[i].set_name(base_name_ + "_right_Mz");
      //face_size = (*region_).ymax() - (*region_).ymin();
      face_size = (*region_).xlen() * (*region_).zlen();

      Jt1_[i].add_dimension("z", (*region_).zlen(), (*gregion).zlen(),
                            (*region_).zstart());
      Jt1_[i].add_dimension("x", (*region_).xlen(), (*gregion).xlen(),
                            (*region_).xstart());

      Jt2_[i].add_dimension("z", (*region_).zlen(), (*gregion).zlen(),
                            (*region_).zstart());
      Jt2_[i].add_dimension("x", (*region_).xlen(), (*gregion).xlen(),
                            (*region_).xstart());

      Mt1_[i].add_dimension("z", (*region_).zlen(), (*gregion).zlen(),
                            (*region_).zstart());
      Mt1_[i].add_dimension("x", (*region_).xlen(), (*gregion).xlen(),
                            (*region_).xstart());

      Mt2_[i].add_dimension("z", (*region_).zlen(), (*gregion).zlen(),
                            (*region_).zstart());
      Mt2_[i].add_dimension("x", (*region_).xlen(), (*gregion).xlen(),
                            (*region_).xstart());
      break;

    case TOP:
      Jt1_[i].set_name(base_name_ + "_top_Jx");
      Jt2_[i].set_name(base_name_ + "_top_Jy");
      Mt1_[i].set_name(base_name_ + "_top_Mx");
      Mt2_[i].set_name(base_name_ + "_top_My");
      //face_size = (*region_).zmax() - (*region_).zmin();
      face_size = (*region_).xlen() * (*region_).ylen();

      Jt1_[i].add_dimension("x", (*region_).xlen(), (*gregion).xlen(),
                            (*region_).xstart());
      Jt1_[i].add_dimension("y", (*region_).ylen(), (*gregion).ylen(),
                            (*region_).ystart());

      Jt2_[i].add_dimension("x", (*region_).xlen(), (*gregion).xlen(),
                            (*region_).xstart());
      Jt2_[i].add_dimension("y", (*region_).ylen(), (*gregion).ylen(),
                            (*region_).ystart());

      Mt1_[i].add_dimension("x", (*region_).xlen(), (*gregion).xlen(),
                            (*region_).xstart());
      Mt1_[i].add_dimension("y", (*region_).ylen(), (*gregion).ylen(),
                            (*region_).ystart());

      Mt2_[i].add_dimension("x", (*region_).xlen(), (*gregion).xlen(),
                            (*region_).xstart());
      Mt2_[i].add_dimension("y", (*region_).ylen(), (*gregion).ylen(),
                            (*region_).ystart());
      break;

    case BOTTOM:
      Jt1_[i].set_name(base_name_ + "_bottom_Jx");
      Jt2_[i].set_name(base_name_ + "_bottom_Jy");
      Mt1_[i].set_name(base_name_ + "_bottom_Mx");
      Mt2_[i].set_name(base_name_ + "_bottom_My");
      //face_size = (*region_).zmax() - (*region_).zmin();
      face_size = (*region_).xlen() * (*region_).ylen();

      Jt1_[i].add_dimension("x", (*region_).xlen(), (*gregion).xlen(),
                            (*region_).xstart());
      Jt1_[i].add_dimension("y", (*region_).ylen(), (*gregion).ylen(),
                            (*region_).ystart());

      Jt2_[i].add_dimension("x", (*region_).xlen(), (*gregion).xlen(),
                            (*region_).xstart());
      Jt2_[i].add_dimension("y", (*region_).ylen(), (*gregion).ylen(),
                            (*region_).ystart());

      Mt1_[i].add_dimension("x", (*region_).xlen(), (*gregion).xlen(),
                            (*region_).xstart());
      Mt1_[i].add_dimension("y", (*region_).ylen(), (*gregion).ylen(),
                            (*region_).ystart());

      Mt2_[i].add_dimension("x", (*region_).xlen(), (*gregion).xlen(),
                            (*region_).xstart());
      Mt2_[i].add_dimension("y", (*region_).ylen(), (*gregion).ylen(),
                            (*region_).ystart());
      break;

    }

    if ((*region_).has_face_data(static_cast<Face>(i)))
    {
      Jt1_data_[i] = new field_t[face_size];
      Jt2_data_[i] = new field_t[face_size];
      Mt1_data_[i] = new field_t[face_size];
      Mt2_data_[i] = new field_t[face_size];

      Jt1_[i].set_ptr(Jt1_data_[i]);
      Jt2_[i].set_ptr(Jt2_data_[i]);
      Mt1_[i].set_ptr(Mt1_data_[i]);
      Mt2_[i].set_ptr(Mt2_data_[i]);

      Jt1_[i].set_num(face_size);
      Jt2_[i].set_num(face_size);
      Mt1_[i].set_num(face_size);
      Mt2_[i].set_num(face_size);

    } else {
      Jt1_[i].set_num(0);
      Jt2_[i].set_num(0);
      Mt1_[i].set_num(0);
      Mt2_[i].set_num(0);
    }
  } // end for (int i = 0; i < 6; i++)

}

void SurfaceCurrentResult::deinit()
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

  if (freqs_)
  {
    delete[] freqs_;
    freqs_ = 0;
  }
}

void SurfaceCurrentResult::calculate_result(const Grid &grid, 
                                            unsigned int time_step)
{

  for (int face_idx = 0; face_idx < 6; face_idx++)
  {
    if (!(*region_).has_face_data(static_cast<Face>(face_idx)))
      continue;

    unsigned int xmin, xmax, ymin, ymax, zmin, zmax;
    
    xmin = (*region_).xmin();
    ymin = (*region_).ymin();
    zmin = (*region_).zmin();
    xmax = (*region_).xmax();
    ymax = (*region_).ymax();
    zmax = (*region_).zmax();

    switch (face_idx)
    {
    case FRONT:
      xmin = xmax - 1;
      calc_currents<YZPlane>(xmin, xmax, ymin, ymax, zmin, zmax, 
                             face_idx, grid);
      break;

    case BACK:
      xmax = xmin + 1;
      calc_currents<YZPlane>(xmin, xmax, ymin, ymax, zmin, zmax, 
                             face_idx, grid);
      break;

    case LEFT:
      ymax = ymin + 1;
      calc_currents<XZPlane>(xmin, xmax, ymin, ymax, zmin, zmax, 
                             face_idx, grid);
      break;

    case RIGHT:
      ymin = ymax - 1;
      calc_currents<XZPlane>(xmin, xmax, ymin, ymax, zmin, zmax, 
                             face_idx, grid);
      break;

    case TOP:
      zmin = zmax - 1;
      calc_currents<XYPlane>(xmin, xmax, ymin, ymax, zmin, zmax, 
                             face_idx, grid);
      break;

    case BOTTOM:
      zmax = zmin + 1;
      calc_currents<XYPlane>(xmin, xmax, ymin, ymax, zmin, zmax, 
                             face_idx, grid);
      break;
    }

  } // end for (int i = 0; i < 6; i++)

}

template<class T>
void SurfaceCurrentResult::calc_currents(unsigned int xmin, 
                                         unsigned int xmax, 
                                         unsigned int ymin, 
                                         unsigned int ymax, 
                                         unsigned int zmin, 
                                         unsigned int zmax, 
                                         int face_idx, 
                                         const Grid &grid)
{
  const field_t *e_t1;
  const field_t *e_t2, *h_t1, *h_t2;

  T p(const_cast<Grid &>(grid)); // EVIL

  field_t *Jt1 = Jt1_data_[face_idx];
  field_t *Jt2 = Jt2_data_[face_idx];

  field_t *Mt1 = Mt1_data_[face_idx];
  field_t *Mt2 = Mt2_data_[face_idx];

  unsigned int idx = 0;

  for (unsigned int i = xmin; i < xmax; i++)
  {
    for (unsigned int j = ymin; j < ymax; j++)
    {
      e_t1 = p.get_e_t1_ptr(i, j, zmin);
      e_t2 = p.get_e_t2_ptr(i, j, zmin);
      h_t1 = p.get_h_t1_ptr(i, j, zmin);
      h_t2 = p.get_h_t2_ptr(i, j, zmin);

      for (unsigned int k = zmin; k < zmax; k++, idx++)
      {
        Jt1[idx] = - *h_t2;
        Jt2[idx] = *h_t1;
        Mt1[idx] = *e_t2;
        Mt2[idx] = - *e_t1;

        h_t2++; h_t1++; e_t1++; e_t2++;
      }
    }
  }
}

ostream& SurfaceCurrentResult::to_string(ostream &os) const
{
  return os << "SurfaceCurrentResult...";
}
