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

#ifndef GRID_UPDATE_H
#define GRID_UPDATE_H

#include "GridUpdateTiling.hh"

/**
 * This is the data class used by the the Electric and Magnetic
 * GridUpdate alogirthms.
 */ 
class GridUpdateData
{
public:
  field_t *ex_;
  field_t *ey_;
  field_t *ez_;
  field_t *hx_;
  field_t *hy_;
  field_t *hz_;

  mat_idx_t *material_;

  mat_coef_t *Ca_;
  mat_coef_t *Cbx_;
  mat_coef_t *Cby_;
  mat_coef_t *Cbz_;
  
  mat_coef_t *Da_;
  mat_coef_t *Dbx_;
  mat_coef_t *Dby_;
  mat_coef_t *Dbz_;

  GridUpdateData(Grid *grid)
    : ex_(grid.get_ex_ptr(0)),
      ey_(grid.get_ey_ptr(0)),
      ez_(grid.get_ez_ptr(0)),
      hx_(grid.get_hx_ptr(0)),
      hy_(grid.get_hy_ptr(0)),
      hz_(grid.get_hz_ptr(0)),
};

/**
 * This class implements an algorithm that updates all three electric
 * field components in once loop, the other three in another. This may
 * improve cache performance over updating each component in it's own
 * loop.
 */ 
class ElectricGridUpdate 
{
public:
  /**
   * Set up the data block as needed... calculate pointers etc. Since
   * pointers aren't really needed, this may go away.
   */ 
  static void pre_z_setup(loop_idx_t x, loop_idx, y
                          loop_idx_t z, Indicies inds, 
                          Data &data)
  {

  }

  /**
   * This algorithm updates the E field components. 
   */
  static void alg(loop_idx_t x, loop_idx, y
                  loop_idx_t z, Indicies inds, 
                  Data &data)
  {
    mat_idx_t mid = data.material_[inds.idx];
    
    data.ex_[inds.idx] = data.Ca_[mid] * data.ex_[inds.idx]
      + data.Cby_[mid] * (data.hz_[inds.idx] - data.hz_[inds.idx_ny])
      + data.Cbz_[mid] * (data.hy_[inds.idx_nz] - data.hy_[inds.idx]);

    data.ey_[inds.idx] = Ca_[mid] * data.ey_[inds.idx]
      + data.Cbz_[mid] * (data.hx_[inds.idx] - data.hx_[inds.idx_ny])
      + data.Cbx_[mid] * (data.hz_[inds.idx_nx] - data.hz_[inds.idx]);    

    data.ez_[inds.idx] = Ca_[mid] * data.ez_[inds.idx]
      + data.Cbx_[mid] * (data.hy_[inds.idx] - data.hy_[inds.idx_nx])
      + data.Cby_[mid] * (data.hx_[inds.idx_ny] - data.hx_[inds.idx]);
  }

};

#endif
