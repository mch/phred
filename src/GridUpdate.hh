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

#ifdef OLD_MATERIAL_DATA
  mat_coef_t * restrict Ca_;
  mat_coef_t * restrict Cbx_;
  mat_coef_t * restrict Cby_;
  mat_coef_t * restrict Cbz_;
  
  mat_coef_t * restrict Da_;
  mat_coef_t * restrict Dbx_;
  mat_coef_t * restrict Dby_;
  mat_coef_t * restrict Dbz_;
#else
  mat_coef_t * restrict C_;
  mat_coef_t * restrict D_;
#endif

  // Need to keep a references so we can call pi()
  Grid &grid_;
  
  GridUpdateData(Grid &grid)
    : ex_(grid.ex_), ey_(grid.ey_), ez_(grid.ez_),
      hx_(grid.hx_), hy_(grid.hy_), hz_(grid.hz_),
      material_(grid.material_), 
#ifdef OLD_MATERIAL_DATA
      Ca_(grid.Ca_),
      Cbx_(grid.Cbx_), Cby_(grid.Cby_), Cbz_(grid.Cbz_),
      Da_(grid.Da_), Dbx_(grid.Dbx_), Dby_(grid.Dby_), Dbz_(grid.Dbz_),
#else
      C_(grid.C_), D_(grid.D_),
#endif    
      grid_(grid)
  {}

#ifdef OLD_MATERIAL_DATA
  /** 
   * Returns the value of the Ca constant for a given material index
   */ 
  inline mat_coef_t get_Ca(mat_idx_t idx)
  { return Ca_[idx]; }

  /**
   * Returns the value of the Cbx constant for a given material index
   */ 
  inline mat_coef_t get_Cbx(mat_idx_t idx)
  { return Cbx_[idx]; }

  /**
   * Returns the value of the Cby constant for a given material index
   */ 
  inline mat_coef_t get_Cby(mat_idx_t idx)
  { return Cby_[idx]; }
  
  /**
   * Returns the value of the Cbz constant for a given material index
   */ 
  inline mat_coef_t get_Cbz(mat_idx_t idx)
  { return Cbz_[idx]; }
  
  /** 
   * Returns the value of the Da constant for a given material index
   */ 
  inline mat_coef_t get_Da(mat_idx_t idx)
  { return Da_[idx]; }

  /**
   * Returns the value of the Dbx constant for a given material index
   */ 
  inline mat_coef_t get_Dbx(mat_idx_t idx)
  { return Dbx_[idx]; }

  /**
   * Returns the value of the Dby constant for a given material index
   */ 
  inline mat_coef_t get_Dby(mat_idx_t idx)
  { return Dby_[idx]; }
  
  /**
   * Returns the value of the Dbz constant for a given material index
   */ 
  inline mat_coef_t get_Dbz(mat_idx_t idx)
  { return Dbz_[idx]; }
#else
  /** 
   * Returns the value of the Ca constant for a given material index
   */ 
  inline mat_coef_t get_Ca(mat_idx_t idx)
  { return C_[idx * 4]; }

  /**
   * Returns the value of the Cbx constant for a given material index
   */ 
  inline mat_coef_t get_Cbx(mat_idx_t idx)
  { return C_[idx * 4 + 1]; }

  /**
   * Returns the value of the Cby constant for a given material index
   */ 
  inline mat_coef_t get_Cby(mat_idx_t idx)
  { return C_[idx * 4 + 2]; }
  
  /**
   * Returns the value of the Cbz constant for a given material index
   */ 
  inline mat_coef_t get_Cbz(mat_idx_t idx)
  { return C_[idx * 4 + 3]; }
  
  /** 
   * Returns the value of the Da constant for a given material index
   */ 
  inline mat_coef_t get_Da(mat_idx_t idx)
  { return D_[idx * 4]; }

  /**
   * Returns the value of the Dbx constant for a given material index
   */ 
  inline mat_coef_t get_Dbx(mat_idx_t idx)
  { return D_[idx * 4 + 1]; }

  /**
   * Returns the value of the Dby constant for a given material index
   */ 
  inline mat_coef_t get_Dby(mat_idx_t idx)
  { return D_[idx * 4 + 2]; }
  
  /**
   * Returns the value of the Dbz constant for a given material index
   */ 
  inline mat_coef_t get_Dbz(mat_idx_t idx)
  { return D_[idx * 4 + 3]; }
#endif
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
    
    data.ex_[pdata.idx] = data.get_Ca(mid) * data.ex_[pdata.idx]
      + data.get_Cby(mid) * (data.hz_[pdata.idx] - data.hz_[pdata.idx_y])
      + data.get_Cbz(mid) * (data.hy_[pdata.idx_z] - data.hy_[pdata.idx]);

    data.ey_[pdata.idx] = data.get_Ca(mid) * data.ey_[pdata.idx]
      + data.get_Cbz(mid) * (data.hx_[pdata.idx] - data.hx_[pdata.idx_z])
      + data.get_Cbx(mid) * (data.hz_[pdata.idx_x] - data.hz_[pdata.idx]);    

    data.ez_[pdata.idx] = data.get_Ca(mid) * data.ez_[pdata.idx]
      + data.get_Cbx(mid) * (data.hy_[pdata.idx] - data.hy_[pdata.idx_x])
      + data.get_Cby(mid) * (data.hx_[pdata.idx_y] - data.hx_[pdata.idx]);

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
    
    data.hx_[pdata.idx] = data.get_Da(mid) * data.hx_[pdata.idx]
      + data.get_Dby(mid) * (data.ez_[pdata.idx] - data.ez_[pdata.idx_y])
      + data.get_Dbz(mid) * (data.ey_[pdata.idx_z] - data.ey_[pdata.idx]);

    data.hy_[pdata.idx] = data.get_Da(mid) * data.hy_[pdata.idx]
      + data.get_Dbz(mid) * (data.ex_[pdata.idx] - data.ex_[pdata.idx_z])
      + data.get_Dbx(mid) * (data.ez_[pdata.idx_x] - data.ez_[pdata.idx]);    

    data.hz_[pdata.idx] = data.get_Da(mid) * data.hz_[pdata.idx]
      + data.get_Dbx(mid) * (data.ey_[pdata.idx] - data.ey_[pdata.idx_x])
      + data.get_Dby(mid) * (data.ex_[pdata.idx_y] - data.ex_[pdata.idx]);

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
    
    data.ex_[pdata.idx] = data.get_Ca(mid) * data.ex_[pdata.idx]
      + data.get_Cby(mid) * (data.hz_[pdata.idx] - data.hz_[pdata.idx_y])
      + data.get_Cbz(mid) * (data.hy_[pdata.idx_z] - data.hy_[pdata.idx]);

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
    
    data.ey_[pdata.idx] = data.get_Ca(mid) * data.ey_[pdata.idx]
      + data.get_Cbz(mid) * (data.hx_[pdata.idx] - data.hx_[pdata.idx_z])
      + data.get_Cbx(mid) * (data.hz_[pdata.idx_x] - data.hz_[pdata.idx]);    

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

    data.ez_[pdata.idx] = data.get_Ca(mid) * data.ez_[pdata.idx]
      + data.get_Cbx(mid) * (data.hy_[pdata.idx] - data.hy_[pdata.idx_x])
      + data.get_Cby(mid) * (data.hx_[pdata.idx_y] - data.hx_[pdata.idx]);    

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
    
    data.hx_[pdata.idx] = data.get_Da(mid) * data.hx_[pdata.idx]
      + data.get_Dby(mid) * (data.ez_[pdata.idx] - data.ez_[pdata.idx_y])
      + data.get_Dbz(mid) * (data.ey_[pdata.idx_z] - data.ey_[pdata.idx]);

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
    
    data.hy_[pdata.idx] = data.get_Da(mid) * data.hy_[pdata.idx]
      + data.get_Dbz(mid) * (data.ex_[pdata.idx] - data.ex_[pdata.idx_z])
      + data.get_Dbx(mid) * (data.ez_[pdata.idx_x] - data.ez_[pdata.idx]);    

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
    

    data.hz_[pdata.idx] = data.get_Da(mid) * data.hz_[pdata.idx]
      + data.get_Dbx(mid) * (data.ey_[pdata.idx] - data.ey_[pdata.idx_x])
      + data.get_Dby(mid) * (data.ex_[pdata.idx_y] - data.ex_[pdata.idx]);

    pdata.idx++;
    pdata.idx_x++;
    pdata.idx_y++;
  }

};


#endif
