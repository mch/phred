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

#ifndef BLOCK_H
#define BLOCK_H

#include "Types.hh"

/**
 * A subset of Yee cells contained by the grid. This object can only
 * be created by a Grid.
 */ 
class Block {
  friend class Grid;
public:

  inline unsigned int xmin() const
  { return xmin_; }
  inline unsigned int ymin() const
  { return ymin_; }
  inline unsigned int zmin() const
  { return zmin_; }

  inline unsigned int xmax() const
  { return xmax_; }
  inline unsigned int ymax() const
  { return ymax_; }
  inline unsigned int zmax() const
  { return zmax_; }

  inline unsigned int xstart() const
  { return start_x_; }

  inline unsigned int ystart() const
  { return start_y_; }

  inline unsigned int zstart() const
  { return start_z_; }

  inline unsigned int xlen() const
  { return len_x_; }

  inline unsigned int ylen() const
  { return len_y_; }

  inline unsigned int zlen() const
  { return len_z_; }

  inline bool has_data() const
  { return has_data_; }

  inline bool has_face_data(Face face) const
  { return faces_[face]; }
  
  inline bool is_global() const
  { return is_global_; }

private:

  Block() // Only Grid can make me!
    : xmin_(0), ymin_(0), zmin_(0), xmax_(0), ymax_(0), zmax_(0),
      start_x_(0), start_y_(0), start_z_(0), len_x_(0), len_y_(0),
      len_z_(0), has_data_(true), is_global_(true)
  {
    for (int i = 0; i < 6; i++)
      faces_[i] = true;
  }

  unsigned int xmin_, ymin_, zmin_;
  unsigned int xmax_, ymax_, zmax_;

  // Starting point of this block in the global grid
  unsigned int start_x_, start_y_, start_z_;

  // Total number of cells, or length, of the total global block, of
  // which this object represents a subset.
  unsigned int len_x_, len_y_, len_z_;

  bool faces_[6]; // True if the face is in the local grid at all. 
  bool has_data_; // True if this block has any data at all in the
                  // local domain.

  bool is_global_; // True if this block of cells is in the global domain. 
};

std::ostream &operator<<(std::ostream &os, const Block &b);

#endif // BLOCK_H
