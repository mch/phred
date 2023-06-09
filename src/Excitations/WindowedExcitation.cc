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

#include "WindowedExcitation.hh"

WindowedExcitation::WindowedExcitation(shared_ptr<Signal> sf) 
  : Excitation(sf), xmin_(0), xmax_(0), ymin_(0), ymax_(0),
    zmin_(0), zmax_(0)
{}

WindowedExcitation::~WindowedExcitation()
{}

void WindowedExcitation::init(const Grid &grid)
{
  Excitation::init(grid);

  if (box_.get())
  {
    point centre = box_->get_centre();
    point size = box_->get_size();
  
    xmin_ = centre.x - size.x / 2;
    ymin_ = centre.y - size.y / 2;
    zmin_ = centre.z - size.z / 2;
      
    xmax_ = centre.x + size.x / 2;
    ymax_ = centre.y + size.y / 2;
    zmax_ = centre.z + size.z / 2;

    lxmin_ = xmin_ + grid.get_deltax() * grid.get_lsx_ol();
    lymin_ = ymin_ + grid.get_deltay() * grid.get_lsy_ol();
    lzmin_ = zmin_ + grid.get_deltaz() * grid.get_lsz_ol();

  }

#ifdef DEBUG
      cerr << "WindowedExcitation::init, xmin_ = " << xmin_
           << ", ymin_ = " << ymin_ << ", zmin_ = " << zmin_
           << ", xmax_ = " << xmax_ << ", ymax_ = " << ymax_
           << ", zmax_ = " << zmax_ << endl;
      cerr << "Local region starts at " << lxmin_ << ", " 
           << lymin_ << ", " << lzmin_ << endl;
#endif
}

void WindowedExcitation::excite(Grid &grid, unsigned int time_step, 
                                FieldType type)
{
  // Find out where we fit in this grid (convert to local coordinates)
  if (type != BOTH && type != type_)
    return;

  if (!region_->has_data())
    return;

  field_t e_time = grid.get_deltat() * time_step;
  field_t h_time = grid.get_deltat() * (time_step - 0.5);
  field_t e_sf = sf_->signal_function(e_time - time_offset_);
  field_t h_sf = sf_->signal_function(h_time - time_offset_);

  if (h_time < time_start_ || (e_time > time_stop_ 
                               && time_stop_ > grid.get_deltat()))
    return;

  field_t e_fld[3];
  field_t h_fld[3];
            
  e_fld[0] = e_sf * polarization_[0];
  e_fld[1] = e_sf * polarization_[1];
  e_fld[2] = e_sf * polarization_[2];

  h_fld[0] = h_sf * polarization_[0];
  h_fld[1] = h_sf * polarization_[1];
  h_fld[2] = h_sf * polarization_[2];

  float fx, fy, fz;
  delta_t dx, dy, dz;

  dx = grid.get_deltax();
  dy = grid.get_deltay();
  dz = grid.get_deltaz();

  if (!soft_) 
  {
    fx = lxmin_;
    for(loop_idx_t i = (*region_).xmin(); i <= (*region_).xmax(); i++, 
          fx += dx)
    {
      fy = lymin_;
      for (loop_idx_t j = (*region_).ymin(); j <= (*region_).ymax(); j++,
             fy += dy)
      {
        fz = lzmin_;
        for (loop_idx_t k = (*region_).zmin(); k <= (*region_).zmax(); k++,
               fz += dz)
        {
          if (type_ == E)
          {
            if (polarization_[0] != 0.0) 
              grid.set_ex(i,j,k, e_fld[0] * window(fx+dx*0.5, fy, fz));

            if (polarization_[1] != 0.0) 
              grid.set_ey(i,j,k, e_fld[1] * window(fx, fy+dy*0.5, fz));

            if (polarization_[2] != 0.0) 
              grid.set_ez(i,j,k, e_fld[2] * window(fx, fy, fz+dz*0.5));
          }

          else if (type_ == H)
          {
            if (polarization_[0] != 0.0) 
              grid.set_hx(i,j,k, h_fld[0] * window(fx, fy+dy*0.5, fz+dz*0.5));

            if (polarization_[1] != 0.0) 
              grid.set_hy(i,j,k, h_fld[1] * window(fx+dx*0.5, fy, fz+dz*0.5));

            if (polarization_[2] != 0.0) 
              grid.set_hz(i,j,k, h_fld[2] * window(fx+dx*0.5, fy+dy*0.5, fz));
          }
        }
      }
    }
  } else {
    fx = lxmin_;
    for(loop_idx_t i = (*region_).xmin(); i <= (*region_).xmax(); i++,
          fx += dx)
    {
      fy = lymin_;
      for (loop_idx_t j = (*region_).ymin(); j <= (*region_).ymax(); j++,
             fy += dy)
      {
        fz = lzmin_;
        for (loop_idx_t k = (*region_).zmin(); k <= (*region_).zmax(); k++,
               fz += dz)
        {
          if (type_ == E)
          {
            if (polarization_[0] != 0.0) 
              grid.set_ex(i,j,k, 
                          window(fx+dx*0.5, fy, fz) * e_fld[0] 
                          + grid.get_ex(i,j,k));

            if (polarization_[1] != 0.0) 
              grid.set_ey(i,j,k, 
                          window(fx, fy+dy*0.5, fz) * e_fld[1] 
                          + grid.get_ey(i,j,k));

            if (polarization_[2] != 0.0) 
              grid.set_ez(i,j,k, 
                          window(fx, fy, fz+dz*0.5) * e_fld[2] 
                          + grid.get_ez(i,j,k));
          }

          else if (type_ == H)
          {
            if (polarization_[0] != 0.0) 
              grid.set_hx(i,j,k, 
                          window(fx, fy+dy*0.5, fz+dz*0.5) * h_fld[0] 
                          + grid.get_hx(i,j,k));
            if (polarization_[1] != 0.0) 
              grid.set_hy(i,j,k, 
                          window(fx+dx*0.5, fy, fz+dz*0.5) * h_fld[1] 
                          + grid.get_hy(i,j,k));
            if (polarization_[2] != 0.0) 
              grid.set_hz(i,j,k, 
                          window(fx+dx*0.5, fy+dy*0.5, fz) * h_fld[2] 
                          + grid.get_hz(i,j,k));
          }
        }
      }
    }
  }
  //cout << "---------------------" << endl;
}

