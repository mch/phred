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

#ifndef GRID_RESULT_H
#define GRID_RESULT_H

#include "Result.hh"

/**
 * This result class makes information about the grid available. 
 */ 
class GridResult : public Result 
{
public:
  GridResult();
  ~GridResult();

  /**
   * Called to perform any initialization that may be required.
   * Converts the region to local grid coordinates and constructs the
   * MPI datatype.
   */
  void init(const Grid &grid);

  /**
   * Called to perform any de-initialization that may be required.
   */
  void deinit();

  /**
   * Print a string representation to an ostream.
   */
  virtual ostream& to_string(ostream &os) const;

private:
  Variable material_ids_; /**< Material id at each cell */ 
  Variable deltaxs_; /**< Cell sizes along the x axis */ 
  Variable deltays_; /**< Cell sizes along the y axis */ 
  Variable deltazs_; /**< Cell sizes along the z axis */ 

  MPI_Datatype datatype_;
};

#endif
