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

#include "GridResult.hh"
#include "../Globals.hh"

GridResult::GridResult()
{}

GridResult::~GridResult()
{}

void GridResult::init(const Grid &grid)
{
  variables_["material_ids"] = &material_ids_;
  //variables_["delta_xs"] = &deltaxs_;
  //variables_["delta_ys"] = &deltays_;
  //variables_["delta_zs"] = &deltazs_;

  if (MPI_SIZE != 1)
    throw DataWriterException("GridResult is only available when running on a single node for now.");
  
  MPI_Type_contiguous(grid.get_ldx() * grid.get_ldy() * grid.get_ldz(), 
                      MPI_UNSIGNED, &datatype_);
  MPI_Type_commit(&datatype_);
  
  material_ids_.add_dimension("x", grid.get_ldx(), grid.get_gdx(), 
                              grid.get_lsx_ol());
  material_ids_.add_dimension("y", grid.get_ldy(), grid.get_gdy(), 
                              grid.get_lsy_ol());
  material_ids_.add_dimension("z", grid.get_ldz(), grid.get_gdz(), 
                              grid.get_lsz_ol());
  
  material_ids_.set_name(base_name_ + "_material_ids");
  material_ids_.set_datatype(datatype_);
  material_ids_.set_ptr(const_cast<unsigned int *>(grid.get_material_ptr(grid_point(0,0,0))));
  
  //deltaxs_.add_dimension

  material_ids_.set_num(1);
  deltaxs_.set_num(0);
  deltays_.set_num(0);
  deltazs_.set_num(0);
}

void GridResult::deinit()
{}

map<string, Variable *> &
GridResult::get_result(const Grid &grid, unsigned int time_step)
{
  if (time_step > 1)
  {
    material_ids_.set_num(0);
    deltaxs_.set_num(0);
    deltays_.set_num(0);
    deltazs_.set_num(0);
  }

  return variables_;
}

