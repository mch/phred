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

#ifndef BLOCK_TILING_H
#define BLOCK_TILING_H

#include "../Block.hh"
#include "../Types.hh"

/**
 * Sample Algorithm Implementation to demonstrate how BlockTiling works. 
 *
 * This algorithm would be used like this:
 *
 * BlockTilingDemo::BlockTilingDemoData data;
 * data.grid = grid;
 * BlockTiling<BlockTilingDemo, BlockTilingDemoData>::loop(block, data);
 *
 * Here block presumably comes from calling something like
 * grid.get_cellset(box_)->get_local_block(), and box_ is a CSGBox.
 *
 * All of the functions must exist, since BlockTiling will try to call
 * them, but they can be empty.
 * 
 * \bug Use traits of some kind to determine if methods are available? 
 */ 
class BlockTilingDemo
{
  /**
   * A data structure that contains the data used by the algorithm. 
   */ 
  typedef struct {
    field_t *ex;
    const Grid &grid;
  } BlockTilingDemoData;
  
  /**
   * This function returns true or false depending on whether or not
   * ::alg should be called. This may be used by BlockTiling to
   * determine if an alorithm should be applied to a specific spot in
   * the grid. This will only be used by the version of BlockTiling
   * which supports multiple algorithms in the inner loop. This test
   * is here instead of in ::alg so that it does not slow down the
   * algorithm when only one algorithm is applied to the entire grid.
   */ 
  static inline bool test(const int &x, const int &y, const int &z, 
                          BlockTilingDemoData &data)
  {
    // Not a real test
    return grid.get_material_id(x,y,z) == 3;
  }

  /** 
   * This function is used to set up items in the data structure such
   * as pointers. This is called before the inner z loop starts. Since
   * memory is contiguous along this dimension, it is useful to
   * compute pointers into the grid just before the inner z loop
   * starts.
   */ 
  static inline void pre_z_setup(const int &x, const int &y, const int &z, 
                                 BlockTilingDemoData &data)
  {
    data.ex = grid.get_pointer(grid_point(x, y, z), data.FC_EX);
  }

  /**
   * A sample algorithm. This is actually very slow because it calls
   * grid.pi() in the inner loop.
   */ 
  static inline void alg(const int &x, const int &y, const int &z, 
                         BlockTilingDemoData &data)
  {
    *(data.ex) = 1.0;
    data.ex++;
  }
};

/**
 * Loops over a Block of cells applying an algorithm.
 *
 * \bug Use Boost's preprocessor library to create a set of templates
 * that can take multiple algorithms.
 */ 
template<class Alg, class Data>
class BlockTiling 
{
public:

  /**
   * Apply the algorithm to the data in the grid in the given Block
   */ 
  static inline void loop(Block &block, Data &data)
  {
    
    for (int i = block.xmin(); i <= block.xmax(); i++)
    {
      for (int j = block.ymin(); j <= block.ymax(); j++)
      {
        Alg::pre_z_setup(x, y, block.zmin(), data);

        for (int k = block.zmin(); k <= block.zmax(); k++)
        {
          Alg::alg(x, y, z, data);
        }
      }
    }
  }

};

#endif // BLOCK_TILING_H
