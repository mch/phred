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


/** \class Excitation
 * \brief A base class for excitations. The default implementation
 * applies a source function to every point in a given
 * region. Subclasses may implement windowing functions or in other
 * ways depend on the position in the grid where the excitation is
 * being applied. 
 */

#ifndef EXCITE_H
#define EXCITE_H

#include "config.h"

#include "Types.hh"
#include "Grid.hh"
#include "SourceFunction.hh"
#include "LifeCycle.hh"

class Excitation : public LifeCycle
{
protected:
  /**
   * Region to apply the source to (in global coordinates)
   */
  region_t region_;

  /**
   * Field polarization to excite, x, y, z
   */
  field_t polarization_[3];

  /**
   * Field to excite, E or H
   */
  FieldType type_;

  /**
   * True if this is a soft source (defaults to false). Source
   * sources are added to the field rather than the field being set
   * to them.
   */
  bool soft_;

  SourceFunction *sf_; /**< Source function object to apply */

public:

  /**
   * Constructs a Excitation object which can apply a source
   * function which depends only on time to a region of a
   * grid. Creates a copy of the source function object to use, so you
   * don't have to hold onto a reference.
   *
   * @param sf SourceFunction object to use as an excitation. 
   */
  Excitation(SourceFunction *sf)
    : type_(E), soft_(false),
      sf_(sf)
  {
    region_.xmin = region_.ymin = region_.zmin = 0;
    region_.xmax = region_.ymax = region_.zmax = 0;
    polarization_[0] = 1;
    polarization_[1] = 0;
    polarization_[2] = 0;
  }

  virtual ~Excitation()
  {}

  /**
   * This function is called to apply this excitation to the grid. 
   *
   * @param grid the grid to operate on
   * @param time_step the time step we are on
   * @param type the field type we are allowed to excite. 
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
            switch (type_) 
            {
            case E:
              if (polarization_[0] != 0.0) grid.set_ex(i,j,k, fld[0]);
              if (polarization_[1] != 0.0) grid.set_ey(i,j,k, fld[1]);
              if (polarization_[2] != 0.0) grid.set_ez(i,j,k, fld[2]);
              break;

            case H:
              if (polarization_[0] != 0.0) grid.set_hx(i,j,k, fld[0]);
              if (polarization_[1] != 0.0) grid.set_hy(i,j,k, fld[1]);
              if (polarization_[2] != 0.0) grid.set_hz(i,j,k, fld[2]);
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
            switch (type_) 
            {
            case E:
              if (polarization_[0] != 0.0) grid.set_ex(i,j,k, 
                                                       fld[0] + grid.get_ex(i,j,k));
              if (polarization_[1] != 0.0) grid.set_ey(i,j,k, 
                                                       fld[1] + grid.get_ey(i,j,k));
              if (polarization_[2] != 0.0) grid.set_ez(i,j,k, 
                                                       fld[2] + grid.get_ez(i,j,k));
              break;

            case H:
              if (polarization_[0] != 0.0) grid.set_hx(i,j,k, 
                                                       fld[0] + grid.get_hx(i,j,k));
              if (polarization_[1] != 0.0) grid.set_hy(i,j,k, 
                                                       fld[1] + grid.get_hy(i,j,k));
              if (polarization_[2] != 0.0) grid.set_hz(i,j,k, 
                                                       fld[2] + grid.get_hz(i,j,k));
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

  /**
   * Set the (rectangular) region to apply the excitation to. The
   * start and end coordinates can range from 0 to global grid size -
   * 1. 
   *
   * @param region the region
   */
  void set_region(region_t region)
  {
    region_ = region;
  }

  /**
   * Set the (rectangular) region to apply the excitation to. The
   * start and end coordinates can range from 0 to global grid size -
   * 1. 
   *
   * @param x_start The starting x coordinate
   * @param x_end The ending x coordinate
   * @param y_start The starting y coordinate
   * @param y_end The ending y coordinate
   * @param z_start The starting z coordinate
   * @param z_end The ending z coordinate
   */ 
  void set_region(unsigned int x_start, unsigned int x_end, 
                  unsigned int y_start, unsigned int y_end, 
                  unsigned int z_start, unsigned int z_end)
  {
    region_.xmin = x_start;
    region_.ymin = y_start;
    region_.zmin = z_start;
    
    region_.xmax = x_end;
    region_.ymax = y_end;
    region_.zmax = z_end;
  }


  /**
   * Set the polarization vector
   *
   * @param x the x component
   * @param y the y component
   * @param z the z component
   */
  void set_polarization(field_t x, field_t y, field_t z)
  {
    polarization_[0] = x;
    polarization_[1] = y;
    polarization_[2] = z;
  }

  /**
   * Set the field to excite. Either E or H. 
   *
   * @param t FieldType
   */
  void set_type(FieldType t)
  {
    type_ = t;
  }

  /**
   * Returns the field type.
   */
  inline FieldType get_type()
  {
    return type_;
  }

  /**
   * Tells if this source is "soft" or not. Soft sources are
   * additive. 
   */
  inline bool get_soft()
  {
    return soft_;
  }

  /**
   * Set the softness of the source.
   * @param soft true if this source should be soft.
   */
  inline void set_soft(bool soft)
  {
    soft_ = soft;
  }

  /**
   * Verify that the excitation region is entierly inside the global grid. 
   */
  virtual void init(const Grid &grid)
  {
    
  }

};

#endif // EXCITE_H
