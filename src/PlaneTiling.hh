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

#ifndef PLANE_TILING_H
#define PLANE_TILING_H

#include "../Block.hh"
#include "../Types.hh"
#include "../GridPlane.hh"

// int's are used as loop counters because that what OpenMP and
// Intel's compiler prefer for vectorization.

typedef struct {
  const field_t *et1, *et2;
  const field_t *ht1, *ht2;

  field_t et1_avg;
  field_t et2_avg;

  field_t ht1_avg;
  field_t ht2_avg;

} Fields_t;

/**
 * This class is a demonstration of how to implement an algorithm that
 * is applied to a specific plane in the grid. 
 */ 
class AlgDemo
{
public:
  
  /** 
   * This class holds input and output data for this particular
   * algorithm. 
   */ 
  class Data
  { 
    friend class AlgDemo;
  public:
    Data(int N)
      : index(N)
    {
      real = new field_t[N];
      imag = new field_t[N];
    }

    ~Data()
    {
      delete[] real;
      delete[] imag;
    }

  private:
    int index;

    field_t *real;
    field_t *imag;
  };
  
  /**
   * This function does any required set up of the data object before
   * looping starts, such as resetting index counters or finding
   * pointers to things
   */ 
  static inline void setup(Data &data)
  {}

  /**
   * static inline allows the compiler to inline the code at the
   * location where it is called, saving the call overhead. 
   *
   * The fields structure contains pointers that can be dereferenced
   * to access the E and H tangential and normal field component
   * averages at the given point. 
   */ 
  static inline void alg(const int &x, const int &y, const int &z, 
                         Fields_t &f, Data &data)
  {
    data.real[data.index] = f.et1_avg * 2.0; //data.cosf;
    data.imag[data.index] = f.et1_avg * 2.0; //data.sinf;
    data.index++;
  }
};

// Forward declaration
template<class Alg, class Data> class PlaneTiling;

/**
 * This class actually does the looping for a YZ plane. 
 */ 
template<class Alg, class Data>
class YZPlaneLoop
{
public:
  friend class PlaneTiling<Alg, Data>;

  static inline void loop(const Grid &grid, Data &data,
                          Block &block, int x)
  {
    int zidx = 0;
    int yidx = 0;

    int zmin = block.zmin();
    int zmax = block.zmax();
    int ymin = block.ymin();
    int ymax = block.ymax();

    // Check against the grid to ensure that we have enough room to
    // compute averages. 
    if (zmax >= grid.get_ldz_sd())
      zmax--;

    if (ymax >= grid.get_ldy_sd())
      ymax--;

    if (x <= 0)
      x++;

    YZPlane plane(const_cast<Grid &>(grid));

    Fields_t f;

    for (yidx = ymin; yidx < ymax; yidx++)
    {
      f.et1 = plane.get_e_t1_ptr(x, yidx, zmin);
      f.et2 = plane.get_e_t2_ptr(x, yidx, zmin);

      f.ht1 = plane.get_h_t1_ptr(x, yidx, zmin);
      f.ht2 = plane.get_h_t2_ptr(x, yidx, zmin);

      for (zidx = zmin; zidx < zmax; zidx++)
      {
        // Impact of copying? Will dead code be removed? 
        f.et1_avg = plane.get_avg_e_t1(x, yidx, zidx);
        f.et2_avg = plane.get_avg_e_t2(x, yidx, zidx);

        f.ht1_avg = plane.get_avg_h_t1(x, yidx, zidx);
        f.ht2_avg = plane.get_avg_h_t2(x, yidx, zidx);

        Alg::alg(x, yidx, zidx, f, data);

        f.et1++; f.et2++; f.ht1++; f.ht2++;
      }
    } // end outer
  }
};

/**
 * This class actually does the looping for a XZ plane. 
 */ 
template<class Alg, class Data>
class XZPlaneLoop
{
public:
  friend class PlaneTiling<Alg, Data>;

  static inline void loop(const Grid &grid, Data &data,
                          Block &block, int y)
  {
    int zidx = 0;
    int xidx = 0;

    int zmin = block.zmin();
    int zmax = block.zmax();
    int xmin = block.xmin();
    int xmax = block.xmax();

    // Check against the grid to ensure that we have enough room to
    // compute averages. 
    if (zmax >= grid.get_ldz_sd())
      zmax--;

    if (xmax >= grid.get_ldx_sd())
      xmax--;

    if (y <= 0)
      y++;

    YZPlane plane(const_cast<Grid &>(grid));

    Fields_t f;

    for (xidx = xmin; xidx < xmax; xidx++)
    {
      f.et1 = plane.get_e_t1_ptr(xidx, y, zmin);
      f.et2 = plane.get_e_t2_ptr(xidx, y, zmin);

      f.ht1 = plane.get_h_t1_ptr(xidx, y, zmin);
      f.ht2 = plane.get_h_t2_ptr(xidx, y, zmin);

      for (zidx = zmin; zidx < zmax; zidx++)
      {
        // Impact of copying? Will dead code be removed? 
        f.et1_avg = plane.get_avg_e_t1(xidx, y, zidx);
        f.et2_avg = plane.get_avg_e_t2(xidx, y, zidx);

        f.ht1_avg = plane.get_avg_h_t1(xidx, y, zidx);
        f.ht2_avg = plane.get_avg_h_t2(xidx, y, zidx);

        Alg::alg(xidx, y, zidx, f, data);

        f.et1++; f.et2++; f.ht1++; f.ht2++;
      }
    } // end outer
  }
};

/**
 * This class actually does the looping for a XY plane. 
 */ 
template<class Alg, class Data>
class XYPlaneLoop
{
public:
  friend class PlaneTiling<Alg, Data>;

  static inline void loop(const Grid &grid, Data &data,
                          Block &block, int z)
  {
    int yidx = 0;
    int xidx = 0;

    int ymin = block.ymin();
    int ymax = block.ymax();
    int xmin = block.xmin();
    int xmax = block.xmax();

    // Check against the grid to ensure that we have enough room to
    // compute averages. 
    if (ymax >= grid.get_ldy_sd())
      ymax--;

    if (xmax >= grid.get_ldx_sd())
      xmax--;

    if (z <= 0)
      z++;

    XYPlane plane(const_cast<Grid &>(grid));

    Fields_t f;

    // The slowest of them all... no contiguous data access possible. 
    for (xidx = xmin; xidx < xmax; xidx++)
    {
      for (yidx = ymin; yidx < ymax; yidx++)
      {
        f.et1 = plane.get_e_t1_ptr(xidx, yidx, z);
        f.et2 = plane.get_e_t2_ptr(xidx, yidx, z);

        f.ht1 = plane.get_h_t1_ptr(xidx, yidx, z);
        f.ht2 = plane.get_h_t2_ptr(xidx, yidx, z);

        // Impact of copying? Will dead code be removed? 
        f.et1_avg = plane.get_avg_e_t1(xidx, yidx, z);
        f.et2_avg = plane.get_avg_e_t2(xidx, yidx, z);

        f.ht1_avg = plane.get_avg_h_t1(xidx, yidx, z);
        f.ht2_avg = plane.get_avg_h_t2(xidx, yidx, z);

        Alg::alg(xidx, yidx, z, f, data);
      }
    } // end outer
  }
};

/**
 * This class performs a tiling operation over a plane. The operation
 * to perform in the inner loop is a function of the class Alg, which
 * much have static inline member functions Alg::setup and Alg::alg,
 * both of which must take a reference to a data object. Alg::alg must
 * have the following prototype:
 *
 * template<class GP, class Data>
 * class Alg {
 * public:
 *   static inline void Alg::alg(Fields_t &f, Data &d) {...}
 * };
 *
 * The Fields_t object holds a set of pointers to the tangential
 * components of the E and H fiels, and a set of 
 *
 * The Data object is just a struct containg any data that the Alg
 * object will need to access. It also serves as the way to get input
 * and output to the Alg::alg function. 
 */ 
template<class Alg, class Data>
class PlaneTiling 
{
public:

  /**
   * Apply the algorithm to the data in the grid in the given range at
   * the face. 
   */ 
  static inline void loop(const Grid &grid, Block &block, 
                          Face &face, Data &data)
  {
    // Get a GridPlane object and set up the loop ranges

    switch (face)
    {
    case FRONT:
      YZPlaneLoop<Alg, Data>::loop(grid, data, block, 
                                   block.xmax() - 1);
      break;

    case BACK:
      YZPlaneLoop<Alg, Data>::loop(grid, data, block, 
                                   block.xmin());
      break;

    case LEFT:
      XZPlaneLoop<Alg, Data>::loop(grid, data, block, 
                                   block.ymin());
      break;

    case RIGHT:
      XZPlaneLoop<Alg, Data>::loop(grid, data, block, 
                                   block.ymax() - 1);      
      break;

    case TOP:
      XYPlaneLoop<Alg, Data>::loop(grid, data, block, 
                                   block.zmax() - 1);
      break;

    case BOTTOM:
      XYPlaneLoop<Alg, Data>::loop(grid, data, block, 
                                   block.zmin());
      break;
    }

  }
};


#endif
