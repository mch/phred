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

/**
 * An excitation that is windowed in 3 space by some windowing function. 
 */
class WindowedExcitation : public Excitation
{
private:
protected:
  /**
   * Subclasses must override this to provide the window
   * function. This is a bit on the wasty side, because it
   * unnecessarilly calculates the x and y parts of the window at
   * every z step. If profiling shows there is too much time wasted
   * here, then optimize this by spliting it into three functions. 
   */
  virtual field_t window(region_t r, unsigned int x, unsigned int y, 
                         unsigned int z) = 0;

public:
  WindowedExcitation(SourceFunction *sf) 
    : Excitation(sf)
  {}

  virtual ~WindowedExcitation()
  {}

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

    region_t r = grid.global_to_local(region_);

    field_t sf = sf_->source_function(grid, time_step);
    field_t fld[3];

    // Bartlett window coefficients
    field_t w = 1;

    fld[0] = sf * polarization_[0];
    fld[1] = sf * polarization_[1];
    fld[2] = sf * polarization_[2];

    if (!soft_) 
    {
      for(unsigned int i = r.xmin; i < r.xmax; i++)
      {
        for (unsigned int j = r.ymin; j < r.ymax; j++)
        {
          for (unsigned int k = r.zmin; k < r.zmax; k++)
          {
            w = window(r, i, j, k);

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
      for(unsigned int i = r.xmin; i < r.xmax; i++)
      {
        for (unsigned int j = r.ymin; j < r.ymax; j++)
        {
          for (unsigned int k = r.zmin; k < r.zmax; k++)
          {
            w = window(r, i, j, k);

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
  }
  
};

#endif // WINDOWED_EXCITATION_H
