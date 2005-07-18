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

#ifndef BLOCK_H
#define BLOCK_H

#include "Types.hh"

/**
 * A subset of Yee cells contained by the grid. This object can only
 * be created by a Grid. 
 *
 * This is only intended to be used by Results and other objects which
 * need to know about the Grid when they have a CSGBox as a starting
 * point. All other uses should use the more mutable region_t defined
 * in Types.hh.
 */ 
class Block {
  friend class Grid;
  friend class CellSet;
public:

  /**
   * First cell in the grid along X which is part of this Block.
   */ 
  inline int xmin() const
  { return xmin_; }

  /**
   * First cell in the grid along Y which is part of this Block.
   */ 
  inline int ymin() const
  { return ymin_; }

  /**
   * First cell in the grid along Z which is part of this Block.
   */ 
  inline int zmin() const
  { return zmin_; }


  /*** WARNING: The ?max() functions below actually return max +
       1. This has to be changed! */ 

  /**
   * The index of the last cell in the grid along X which is part of
   * this Block.
   */ 
  inline int xmax() const
  { return xmax_; }

  /**
   * The index of the last cell in the grid along Y which is part of
   * this Block.
   */ 
  inline int ymax() const
  { return ymax_; }

  /**
   * The index of the last cell in the grid along Z which is part of
   * this Block.
   */ 
  inline int zmax() const
  { return zmax_; }

  /**
   * Returns the offset of this local grid Block from the start of the
   * global Block. This is used to decide where in the global block
   * this local contribution belongs. Usually for results. 
   */ 
  inline int xoffset() const
  { return xoffset_; }

  inline int yoffset() const
  { return yoffset_; }

  inline int zoffset() const
  { return zoffset_; }

  /** 
   * The length of this Block along each axis. 
   */ 
  inline int xlen() const
  { return xlen_; }

  inline int ylen() const
  { return ylen_; }

  inline int zlen() const
  { return zlen_; }

  /**
   * True if this Block has any data. If a CSGObject is entirely
   * outside the grid and is_global() == true, or if the CSGObject
   * occupies zero cells on this rank and is_globa() == false, then
   * this will return false.
   */ 
  inline bool has_data() const
  { return has_data_; }

  /**
   * Returns true if the given face of this Block is within the local
   * grid.
   */
  inline bool has_face_data(Face face) const
  { return faces_[face]; }
  
  /**
   * Returns true if this Block represents a set of cells in the
   * global computational domain.
   */ 
  inline bool is_global() const
  { return is_global_; }

  /**
   * Returns the number of cells this Block occupies.
   */ 
  inline int volume() const
  { return (xmax_ - xmin_) * (ymax_ - ymin_) * (zmax_ - zmin_); }

private:

  Block() // Only Grid can make me!
    : xmin_(0), ymin_(0), zmin_(0), xmax_(0), ymax_(0), zmax_(0),
      xoffset_(0), yoffset_(0), zoffset_(0),
      xlen_(0), ylen_(0), zlen_(0), 
      has_data_(true), is_global_(true)
  {
    for (int i = 0; i < 6; i++)
      faces_[i] = true;
  }

  int xmin_, ymin_, zmin_;
  int xmax_, ymax_, zmax_;

  // Offset of the start of this local block from the start of the
  // global block.
  int xoffset_, yoffset_, zoffset_;
  
  // Total number of cells, or length, of the total global block, of
  // which this object represents a subset.
  int xlen_, ylen_, zlen_;

  bool faces_[6]; // True if the face is in the local grid at all. 
  bool has_data_; // True if this block has any data at all in the
                  // local domain.

  bool is_global_; // True if this block of cells is in the global domain. 
};

std::ostream &operator<<(std::ostream &os, const Block &b);

#endif // BLOCK_H
