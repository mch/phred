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

#ifndef WINDOWED_EXCITATION_H
#define WINDOWED_EXCITATION_H

#include "Excitation.hh"
#include "../Globals.hh"

/**
 * An excitation that is windowed in 3 space by some windowing function. 
 */
class WindowedExcitation : public Excitation
{
private:
protected:

  // Box ranges in real global coordinates. For help in calculating
  // the window.
  float xmin_, xmax_, ymin_, ymax_, zmin_, zmax_;

  // The part of the box that is in the local domain. For looping
  // across the cells contained in the box.
  float lxmin_, lymin_, lzmin_;

public:
  WindowedExcitation(shared_ptr<SourceFunction> sf) 
    : Excitation(sf), xmin_(0), xmax_(0), ymin_(0), ymax_(0),
      zmin_(0), zmax_(0)
  {}

  virtual ~WindowedExcitation()
  {}

  /**
   * init the windowing function.
   */ 
  virtual void init(const Grid &grid)
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

//     cerr << "WindowedExcitation::init, xmin_ = " << xmin_
//          << ", ymin_ = " << ymin_ << ", zmin_ = " << zmin_
//          << ", xmax_ = " << xmax_ << ", ymax_ = " << ymax_
//          << ", zmax_ = " << zmax_ << endl;
//     cerr << "Local region starts at " << lxmin_ << ", " 
//          << lymin_ << ", " << lzmin_ << endl;
  }

  /**
   * Subclasses must override this to provide the window
   * function. This is a bit on the wasty side, because it
   * unnecessarilly calculates the x and y parts of the window at
   * every z step. If profiling shows there is too much time wasted
   * here, then optimize this by spliting it into three functions. 
   *
   * @param x the x position in the global real coordinate system
   * @param y the y position in the global real coordinate system
   * @param z the z position in the global real coordinate system
   * @return a value describing the strength of the field at the given point. 
   */
  virtual field_t window(float x, float y, float z) = 0;

  /**
   * This applied the windowed excitation to the grid.
   *
   * @param grid the grid to mess with
   * @param time_step the time step we are on
   * @param type the field components to excite (E or H)
   */
  virtual void excite(Grid &grid, unsigned int time_step, 
                      FieldType type)
  {
    // Find out where we fit in this grid (convert to local coordinates)
    if (type != BOTH && type != type_)
      return;

    field_t e_sf = sf_->source_function(grid, time_step);
    field_t h_sf = sf_->source_function(grid, time_step - 0.5);

    field_t e_fld[3];
    field_t h_fld[3];
            
    // Bartlett window coefficients
    field_t w = 1;

    e_fld[0] = e_sf * polarization_[0];
    e_fld[1] = e_sf * polarization_[1];
    e_fld[2] = e_sf * polarization_[2];

    h_fld[0] = h_sf * polarization_[0];
    h_fld[1] = h_sf * polarization_[1];
    h_fld[2] = h_sf * polarization_[2];

    float fx, fy, fz;

    //cout << "Windowing function on " << MPI_RANK 
    //     << ": ---------------------" << endl;
    if (!soft_) 
    {
      fx = lxmin_;
      for(unsigned int i = (*region_).xmin(); i < (*region_).xmax(); i++, 
            fx += grid.get_deltax())
      {
        fy = lymin_;
        for (unsigned int j = (*region_).ymin(); j < (*region_).ymax(); j++,
               fy += grid.get_deltay())
        {
          fz = lzmin_;
          for (unsigned int k = (*region_).zmin(); k < (*region_).zmax(); k++,
                 fz += grid.get_deltaz())
          {
            w = window(fx, fy, fz);
            //cout << fx << " " << fy << " " << fz << " " << w << endl;

            switch (type_) 
            {
            case E:
              if (polarization_[0] != 0.0) grid.set_ex(i,j,k, e_fld[0] * w);
              if (polarization_[1] != 0.0) grid.set_ey(i,j,k, e_fld[1] * w);
              if (polarization_[2] != 0.0) grid.set_ez(i,j,k, e_fld[2] * w);
              break;

            case H:
              if (polarization_[0] != 0.0) grid.set_hx(i,j,k, h_fld[0] * w);
              if (polarization_[1] != 0.0) grid.set_hy(i,j,k, h_fld[1] * w);
              if (polarization_[2] != 0.0) grid.set_hz(i,j,k, h_fld[2] * w);
              break;

            case BOTH: // Isn't meant for Excitations.
              throw std::exception();
              break; 
            }
          }
        }
      }
    } else {
      fx = lxmin_;
      for(unsigned int i = (*region_).xmin(); i < (*region_).xmax(); i++,
            fx += grid.get_deltax())
      {
        fy = lymin_;
        for (unsigned int j = (*region_).ymin(); j < (*region_).ymax(); j++,
               fy += grid.get_deltay())
        {
          fz = lzmin_;
          for (unsigned int k = (*region_).zmin(); k < (*region_).zmax(); k++,
                 fz += grid.get_deltaz())
          {
            w = window(fx, fy, fz);
            //cout << fx << " " << fy << " " << fz << " " << w << endl;

            switch (type_) 
            {
            case E:
              if (polarization_[0] != 0.0) grid.set_ex(i,j,k, 
                                                       w * e_fld[0] 
                                                       + grid.get_ex(i,j,k));
              if (polarization_[1] != 0.0) grid.set_ey(i,j,k, 
                                                       w * e_fld[1] 
                                                       + grid.get_ey(i,j,k));
              if (polarization_[2] != 0.0) grid.set_ez(i,j,k, 
                                                       w * e_fld[2] 
                                                       + grid.get_ez(i,j,k));
              break;

            case H:
              if (polarization_[0] != 0.0) grid.set_hx(i,j,k, 
                                                       w * h_fld[0] 
                                                       + grid.get_hx(i,j,k));
              if (polarization_[1] != 0.0) grid.set_hy(i,j,k, 
                                                       w * h_fld[1] 
                                                       + grid.get_hy(i,j,k));
              if (polarization_[2] != 0.0) grid.set_hz(i,j,k, 
                                                       w * h_fld[2] 
                                                       + grid.get_hz(i,j,k));
              break;

            case BOTH: // Isn't meant for Excitations.
              throw std::exception();
              break; 
            }
          }
        }
      }
    }
    //cout << "---------------------" << endl;
  }
  
};

#endif // WINDOWED_EXCITATION_H
