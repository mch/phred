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

  // Box ranges in real global coordinates. 
  float xmin_, xmax_, ymin_, ymax_, zmin_, zmax_;

public:
  WindowedExcitation(SourceFunction *sf) 
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
    }
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

    field_t sf = sf_->source_function(grid, time_step);
    field_t fld[3];
            
    // Bartlett window coefficients
    field_t w = 1;

    fld[0] = sf * polarization_[0];
    fld[1] = sf * polarization_[1];
    fld[2] = sf * polarization_[2];

    float fx, fy, fz;

    //cout << "Windowing function on " << MPI_RANK 
    //     << ": ---------------------" << endl;
    if (!soft_) 
    {
      fx = grid.get_deltax() * (grid.get_lsx_ol() + region_.xmin);
      for(unsigned int i = region_.xmin; i < region_.xmax; i++, 
            fx += grid.get_deltax())
      {
        fy = grid.get_deltay() * (grid.get_lsy_ol() + region_.ymin);
        for (unsigned int j = region_.ymin; j < region_.ymax; j++,
               fy += grid.get_deltay())
        {
          fz = grid.get_deltaz() * (grid.get_lsz_ol() + region_.zmin);
          for (unsigned int k = region_.zmin; k < region_.zmax; k++,
                 fz += grid.get_deltaz())
          {
            w = window(fx, fy, fz);
            //cout << gi << " " << gj << " " << gk << " " << w << endl;

            switch (type_) 
            {
            case E:
              if (polarization_[0] != 0.0) grid.set_ex(i,j,k, fld[0] * w);
              if (polarization_[1] != 0.0) grid.set_ey(i,j,k, fld[1] * w);
              if (polarization_[2] != 0.0) grid.set_ez(i,j,k, fld[2] * w);
              break;

            case H:
              if (polarization_[0] != 0.0) grid.set_hx(i,j,k, fld[0] * w);
              if (polarization_[1] != 0.0) grid.set_hy(i,j,k, fld[1] * w);
              if (polarization_[2] != 0.0) grid.set_hz(i,j,k, fld[2] * w);
              break;

            case BOTH: // Isn't meant for Excitations.
              throw std::exception();
              break; 
            }
          }
        }
      }
    } else {
      fx = grid.get_deltax() * (grid.get_lsx_ol() + region_.xmin);
      for(unsigned int i = region_.xmin; i < region_.xmax; i++,
            fx += grid.get_deltax())
      {
        fy = grid.get_deltay() * (grid.get_lsy_ol() + region_.ymin);
        for (unsigned int j = region_.ymin; j < region_.ymax; j++,
               fy += grid.get_deltay())
        {
          fz = grid.get_deltaz() * (grid.get_lsz_ol() + region_.zmin);
          for (unsigned int k = region_.zmin; k < region_.zmax; k++,
                 fz += grid.get_deltaz())
          {
            w = window(fx, fy, fz);
            //cout << gi << " " << gj << " " << gk << " " << w << endl;

            switch (type_) 
            {
            case E:
              if (polarization_[0] != 0.0) grid.set_ex(i,j,k, 
                                                       w * fld[0] 
                                                       + grid.get_ex(i,j,k));
              if (polarization_[1] != 0.0) grid.set_ey(i,j,k, 
                                                       w * fld[1] 
                                                       + grid.get_ey(i,j,k));
              if (polarization_[2] != 0.0) grid.set_ez(i,j,k, 
                                                       w * fld[2] 
                                                       + grid.get_ez(i,j,k));
              break;

            case H:
              if (polarization_[0] != 0.0) grid.set_hx(i,j,k, 
                                                       w * fld[0] 
                                                       + grid.get_hx(i,j,k));
              if (polarization_[1] != 0.0) grid.set_hy(i,j,k, 
                                                       w * fld[1] 
                                                       + grid.get_hy(i,j,k));
              if (polarization_[2] != 0.0) grid.set_hz(i,j,k, 
                                                       w * fld[2] 
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
