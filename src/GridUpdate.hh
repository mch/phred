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
  // If we are keeping a reference, it is really necessary to copy all
  // these pointers? By adding the restrict keyword, we tell the
  // compiler that none of these arrays overlap, and that more
  // aggressive optimizations (and SIMD parallism) may be possible...
  field_t * restrict ex_;
  field_t * restrict ey_;
  field_t * restrict ez_;
  field_t * restrict hx_;
  field_t * restrict hy_;
  field_t * restrict hz_;

  mat_idx_t * restrict material_;

  mat_coef_t * restrict Ca_;
  mat_coef_t * restrict Cbx_;
  mat_coef_t * restrict Cby_;
  mat_coef_t * restrict Cbz_;
  
  mat_coef_t * restrict Da_;
  mat_coef_t * restrict Dbx_;
  mat_coef_t * restrict Dby_;
  mat_coef_t * restrict Dbz_;

  // Need to keep a references so we can call pi()
  Grid &grid_;
  
  GridUpdateData(Grid &grid)
    : ex_(grid.ex_), ey_(grid.ey_), ez_(grid.ez_),
      hx_(grid.hx_), hy_(grid.hy_), hz_(grid.hz_),
      material_(grid.material_), Ca_(grid.Ca_),
      Cbx_(grid.Cbx_), Cby_(grid.Cby_), Cbz_(grid.Cbz_),
      Da_(grid.Da_), Dbx_(grid.Dbx_), Dby_(grid.Dby_), Dbz_(grid.Dbz_),
      grid_(grid)
  {}
};

/**
 * This class holds data that must be private to OpenMP threads,
 * such as indicies into arrays.
 * 
 * First should be at (i,j,k), last three are plus or minus 1 along
 * each axis, depending on the field being updated.
 */ 
class PrivateGridUpdateData
{
public:
  int idx, idx_x, idx_y, idx_z;
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
   * Set up the data block as needed... 
   */ 
  static void pre_z_setup(const loop_idx_t &i, const loop_idx_t &j,
                          const loop_idx_t &k, GridUpdateData &data,
                          PrivateGridUpdateData &pdata)
  { 
    pdata.idx = data.grid_.pi(i, j, k);

    pdata.idx_x = data.grid_.pi(i-1, j, k);
    pdata.idx_y = data.grid_.pi(i, j-1, k);
    pdata.idx_z = pdata.idx - 1;
  }

  /**
   * This algorithm updates the E field components. 
   */
  static void alg(const loop_idx_t &x, const loop_idx_t &y,
                  const loop_idx_t &z, GridUpdateData &data,
                  PrivateGridUpdateData &pdata)
  {
    mat_idx_t mid = data.material_[pdata.idx];
    
    data.ex_[pdata.idx] = data.Ca_[mid] * data.ex_[pdata.idx]
      + data.Cby_[mid] * (data.hz_[pdata.idx] - data.hz_[pdata.idx_y])
      + data.Cbz_[mid] * (data.hy_[pdata.idx_z] - data.hy_[pdata.idx]);

    data.ey_[pdata.idx] = data.Ca_[mid] * data.ey_[pdata.idx]
      + data.Cbz_[mid] * (data.hx_[pdata.idx] - data.hx_[pdata.idx_z])
      + data.Cbx_[mid] * (data.hz_[pdata.idx_x] - data.hz_[pdata.idx]);    

    data.ez_[pdata.idx] = data.Ca_[mid] * data.ez_[pdata.idx]
      + data.Cbx_[mid] * (data.hy_[pdata.idx] - data.hy_[pdata.idx_x])
      + data.Cby_[mid] * (data.hx_[pdata.idx_y] - data.hx_[pdata.idx]);

    pdata.idx++;
    pdata.idx_x++;
    pdata.idx_y++;
    pdata.idx_z++;
  }

};

/**
 * This class implements an algorithm that updates all three magnetic
 * field components in one loop, the other three in another. This may
 * improve cache performance over updating each component in it's own
 * loop.
 */ 
class MagneticGridUpdate 
{
public:
  /**
   * Set up the data block as needed... 
   */ 
  static void pre_z_setup(const loop_idx_t &i, const loop_idx_t &j,
                          const loop_idx_t &k, GridUpdateData &data,
                          PrivateGridUpdateData &pdata)
  { 
    pdata.idx = data.grid_.pi(i, j, k);

    pdata.idx_x = data.grid_.pi(i+1, j, k);
    pdata.idx_y = data.grid_.pi(i, j+1, k);
    pdata.idx_z = pdata.idx + 1;
  }

  /**
   * This algorithm updates the H field components. 
   */
  static void alg(const loop_idx_t &x, const loop_idx_t &y,
                  const loop_idx_t &z, GridUpdateData &data,
                  PrivateGridUpdateData &pdata)
  {
    mat_idx_t mid = data.material_[pdata.idx];
    
    data.hx_[pdata.idx] = data.Da_[mid] * data.hx_[pdata.idx]
      + data.Dby_[mid] * (data.ez_[pdata.idx] - data.ez_[pdata.idx_y])
      + data.Dbz_[mid] * (data.ey_[pdata.idx_z] - data.ey_[pdata.idx]);

    data.hy_[pdata.idx] = data.Da_[mid] * data.hy_[pdata.idx]
      + data.Dbz_[mid] * (data.ex_[pdata.idx] - data.ex_[pdata.idx_z])
      + data.Dbx_[mid] * (data.ez_[pdata.idx_x] - data.ez_[pdata.idx]);    

    data.hz_[pdata.idx] = data.Da_[mid] * data.hz_[pdata.idx]
      + data.Dbx_[mid] * (data.ey_[pdata.idx] - data.ey_[pdata.idx_x])
      + data.Dby_[mid] * (data.ex_[pdata.idx_y] - data.ex_[pdata.idx]);

    pdata.idx++;
    pdata.idx_x++;
    pdata.idx_y++;
    pdata.idx_z++;
  }

};

/**
 * This class updates on the Ex field only. 
 */ 
class ExGridUpdate 
{
public:
  /**
   * Set up the data block as needed... 
   */ 
  static void pre_z_setup(const loop_idx_t &i, const loop_idx_t &j,
                          const loop_idx_t &k, GridUpdateData &data,
                          PrivateGridUpdateData &pdata)
  { 
    pdata.idx = data.grid_.pi(i, j, k);

    pdata.idx_y = data.grid_.pi(i, j-1, k);
    pdata.idx_z = pdata.idx - 1;
  }

  /**
   * This algorithm updates the Ex field. 
   */
  static void alg(const loop_idx_t &x, const loop_idx_t &y,
                  const loop_idx_t &z, GridUpdateData &data,
                  PrivateGridUpdateData &pdata)
  {
    mat_idx_t mid = data.material_[pdata.idx];
    
    data.ex_[pdata.idx] = data.Ca_[mid] * data.ex_[pdata.idx]
      + data.Cby_[mid] * (data.hz_[pdata.idx] - data.hz_[pdata.idx_y])
      + data.Cbz_[mid] * (data.hy_[pdata.idx_z] - data.hy_[pdata.idx]);

    pdata.idx++;
    pdata.idx_y++;
    pdata.idx_z++;
  }

};

/**
 * This class updates on the Ey field only. 
 */ 
class EyGridUpdate 
{
public:
  /**
   * Set up the data block as needed... 
   */ 
  static void pre_z_setup(const loop_idx_t &i, const loop_idx_t &j,
                          const loop_idx_t &k, GridUpdateData &data,
                          PrivateGridUpdateData &pdata)
  { 
    pdata.idx = data.grid_.pi(i, j, k);

    pdata.idx_x = data.grid_.pi(i-1, j, k);
    pdata.idx_z = pdata.idx - 1;
  }

  /**
   * This algorithm updates the Ey field. 
   */
  static void alg(const loop_idx_t &x, const loop_idx_t &y,
                  const loop_idx_t &z, GridUpdateData &data,
                  PrivateGridUpdateData &pdata)
  {
    mat_idx_t mid = data.material_[pdata.idx];
    
    data.ey_[pdata.idx] = data.Ca_[mid] * data.ey_[pdata.idx]
      + data.Cbz_[mid] * (data.hx_[pdata.idx] - data.hx_[pdata.idx_z])
      + data.Cbx_[mid] * (data.hz_[pdata.idx_x] - data.hz_[pdata.idx]);    

    pdata.idx++;
    pdata.idx_x++;
    pdata.idx_z++;
  }

};

/**
 * This class updates on the Ez field only. 
 */ 
class EzGridUpdate 
{
public:
  /**
   * Set up the data block as needed... 
   */ 
  static void pre_z_setup(const loop_idx_t &i, const loop_idx_t &j,
                          const loop_idx_t &k, GridUpdateData &data,
                          PrivateGridUpdateData &pdata)
  { 
    pdata.idx = data.grid_.pi(i, j, k);

    pdata.idx_x = data.grid_.pi(i-1, j, k);
    pdata.idx_y = data.grid_.pi(i, j-1, k);
  }

  /**
   * This algorithm updates the Ez field. 
   */
  static void alg(const loop_idx_t &x, const loop_idx_t &y,
                  const loop_idx_t &z, GridUpdateData &data,
                  PrivateGridUpdateData &pdata)
  {
    mat_idx_t mid = data.material_[pdata.idx];

    data.ez_[pdata.idx] = data.Ca_[mid] * data.ez_[pdata.idx]
      + data.Cbx_[mid] * (data.hy_[pdata.idx] - data.hy_[pdata.idx_x])
      + data.Cby_[mid] * (data.hx_[pdata.idx_y] - data.hx_[pdata.idx]);    

    pdata.idx++;
    pdata.idx_x++;
    pdata.idx_y++;
  }

};

/**
 * This class updates on the Hx field only. 
 */ 
class HxGridUpdate 
{
public:
  /**
   * Set up the data block as needed... 
   */ 
  static void pre_z_setup(const loop_idx_t &i, const loop_idx_t &j,
                          const loop_idx_t &k, GridUpdateData &data,
                          PrivateGridUpdateData &pdata)
  { 
    pdata.idx = data.grid_.pi(i, j, k);

    pdata.idx_y = data.grid_.pi(i, j+1, k);
    pdata.idx_z = pdata.idx + 1;
  }

  /**
   * This algorithm updates the Hx field. 
   */
  static void alg(const loop_idx_t &x, const loop_idx_t &y,
                  const loop_idx_t &z, GridUpdateData &data,
                  PrivateGridUpdateData &pdata)
  {
    mat_idx_t mid = data.material_[pdata.idx];
    
    data.hx_[pdata.idx] = data.Da_[mid] * data.hx_[pdata.idx]
      + data.Dby_[mid] * (data.ez_[pdata.idx] - data.ez_[pdata.idx_y])
      + data.Dbz_[mid] * (data.ey_[pdata.idx_z] - data.ey_[pdata.idx]);

    pdata.idx++;
    pdata.idx_y++;
    pdata.idx_z++;
  }

};

/**
 * This class updates on the Hy field only. 
 */ 
class HyGridUpdate 
{
public:
  /**
   * Set up the data block as needed... 
   */ 
  static void pre_z_setup(const loop_idx_t &i, const loop_idx_t &j,
                          const loop_idx_t &k, GridUpdateData &data,
                          PrivateGridUpdateData &pdata)
  { 
    pdata.idx = data.grid_.pi(i, j, k);

    pdata.idx_x = data.grid_.pi(i+1, j, k);
    pdata.idx_z = pdata.idx + 1;
  }

  /**
   * This algorithm updates the Hy field. 
   */
  static void alg(const loop_idx_t &x, const loop_idx_t &y,
                  const loop_idx_t &z, GridUpdateData &data,
                  PrivateGridUpdateData &pdata)
  {
    mat_idx_t mid = data.material_[pdata.idx];
    
    data.hy_[pdata.idx] = data.Da_[mid] * data.hy_[pdata.idx]
      + data.Dbz_[mid] * (data.ex_[pdata.idx] - data.ex_[pdata.idx_z])
      + data.Dbx_[mid] * (data.ez_[pdata.idx_x] - data.ez_[pdata.idx]);    

    pdata.idx++;
    pdata.idx_x++;
    pdata.idx_z++;
  }

};

/**
 * This class updates on the Hz field only. 
 */ 
class HzGridUpdate 
{
public:
  /**
   * Set up the data block as needed... 
   */ 
  static void pre_z_setup(const loop_idx_t &i, const loop_idx_t &j,
                          const loop_idx_t &k, GridUpdateData &data,
                          PrivateGridUpdateData &pdata)
  { 
    pdata.idx = data.grid_.pi(i, j, k);

    pdata.idx_x = data.grid_.pi(i+1, j, k);
    pdata.idx_y = data.grid_.pi(i, j+1, k);
  }

  /**
   * This algorithm updates the Hz field. 
   */
  static void alg(const loop_idx_t &x, const loop_idx_t &y,
                  const loop_idx_t &z, GridUpdateData &data,
                  PrivateGridUpdateData &pdata)
  {
    mat_idx_t mid = data.material_[pdata.idx];
    

    data.hz_[pdata.idx] = data.Da_[mid] * data.hz_[pdata.idx]
      + data.Dbx_[mid] * (data.ey_[pdata.idx] - data.ey_[pdata.idx_x])
      + data.Dby_[mid] * (data.ex_[pdata.idx_y] - data.ex_[pdata.idx]);

    pdata.idx++;
    pdata.idx_x++;
    pdata.idx_y++;
  }

};


#endif
