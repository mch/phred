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

#ifndef GRID_UPDATE_TILING_H
#define GRID_UPDATE_TILING_H

#include "../Block.hh"
#include "../Types.hh"
#include "../Grid.hh"

/**
 * This class provides tiling methods which loop over a Block of Grid
 * cells for the purpose of applying an FDTD update equation. Each
 * method is specialized for a specific field component, so the
 * specific ranges for the individual field components are
 * encapsulated in the loop ranges. The loops are also specialized in
 * that the required field components are set up and passed to the
 * algorithm.
 *
 * \bug Use Boost's preprocessor library to create a set of templates
 * that can take multiple algorithms.
 */
class GridUpdateTiling {
public:
  
  /**
   * Ensure the given Block will not cause the algorithm to reach
   * outside of the memory allocated to the Grid.
   */
  static inline void check_block(Block &block)
  {
    
  }

  /**
   * Apply an algorithm to the data in the grid in the given
   * Block. The given algorithm MUST implement an update equation for
   * the Ex field component.
   */ 
  static inline void ex_loop(Block &block, Data &data)
  {
    check_block(block);

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

#endif // GRID_UPDATE_TILING_H
