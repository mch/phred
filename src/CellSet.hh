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

#ifndef CELL_SET_H
#define CELL_SET_H

#include "Block.hh"

/**
 * This class represents a set of cells in the computational
 * domain. It includes information about the set of cells in the
 * global computational domain, the local domain WITH ghost cells, and
 * the local domain WITHOUT ghost cells. 
 * 
 * The position in the global computational domain is calculated by
 * the Grid from some CSGObject. The Grid then calculates the local
 * Blocks from the global Block by taking the sub-domaining into
 * account. 
 */ 
class CellSet {
public:
  /**
   * Returns the Block descibing this cell set with respect to the
   * global computational domain. 
   */
  inline shared_ptr<Block> get_global_block() const
  { return global_; }

  /**
   * Returns the Block describing this cell set with respect to the
   * local computational domain. Ghost cells are NOT included. 
   */
  inline shared_ptr<Block> get_local_block() const
  { return local_; }

  /**
   * Returns the Block describing this cell set with respect to the
   * local computational domain. Ghost cells ARE included. 
   */
  shared_ptr<Block> get_local_block_with_ghost() const
  { return local_ghost_; }

  ~CellSet() {}

private:
  friend class Grid;
  CellSet()
    : global_(shared_ptr<Block>(new Block)),
      local_(shared_ptr<Block>(new Block)),
      local_ghost_(shared_ptr<Block>(new Block))
  {}

  /** 
   * This cellset in the global computational domain. 
   */ 
  shared_ptr<Block> global_;

  /**
   * This cellset in the local computational domain, excluding ghost
   * cells. 
   */ 
  shared_ptr<Block> local_;

  /**
   * This cellset in the local computational domain, including ghost
   * cells. 
   */ 
  shared_ptr<Block> local_ghost_;

};

#endif // CELL_SET_H
