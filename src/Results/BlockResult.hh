/* 
   phred - Phred is a parallel finite difference time domain
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

#ifndef BLOCK_RESULT_H
#define BLOCK_RESULT_H

#include "Result.hh"
#include "Types.hh"

#include <mpi.h>

/**
 * This result produces a block of data within the grid. Generally it
 * will only return data for only one field component, but it can
 * optionally return them all. 
 */
class BlockResult : public Result
{
private:
protected:
  region_t region_;
  FieldComponent field_comp_;

  bool init_; /**< Set to true after init() has been called. */

  Variable var_; /**< Our variable */

  MPI_Datatype datatype_;

public:
  BlockResult();
  BlockResult(region_t r, FieldComponent field_comp = FC_EY);
  ~BlockResult();

  /**
   * Return an MPI datatype with a pointer to the block of data. 
   *
   * @param grid a reference to a Grid object
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  map<string, Variable *> &get_result(const Grid &grid, 
                                      unsigned int time_step);

  /**
   * Set the region to return. 
   *
   * @param region
   */
  inline void set_region(region_t region)
  {
    region_ = region;
  }

  /**
   * Called to perform any initialization that may be required.
   * Converts the region to local grid coordinates and constructs the
   * MPI datatype.
   */
  virtual void init(const Grid &grid);

  /**
   * Called to perform any de-initialization that may be required.
   */
  virtual void deinit();

  /**
   * Returns the region this result deals with. 
   *
   * @return the region
   */
  inline region_t get_region()
  {
    return region_;
  }

  /**
   * Set the field compoment to operate on
   *
   * @param field_comp field component
   */
  inline void set_field(FieldComponent field_comp)
  {
    field_comp_ = field_comp;
  }

  /**
   * Returns the field component we are returning
   *
   * @return a field component name
   */
  inline FieldComponent get_field()
  {
    return field_comp_;
  }

};

#endif // BLOCK_RESULT_H
