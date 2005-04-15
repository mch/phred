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

#include "PropertiesResult.hh"

void PropertiesResult::init(const Grid &grid)
{
  dx_var_.set_name(base_name_ + "_dx");
  dy_var_.set_name(base_name_ + "_dy");
  dz_var_.set_name(base_name_ + "_dz");
  dt_var_.set_name(base_name_ + "_dt");
  ts_var_.set_name(base_name_ + "_ts");

  dx_var_.add_dimension("x", 1, 1, 0);
  dy_var_.add_dimension("x", 1, 1, 0);
  dz_var_.add_dimension("x", 1, 1, 0);
  dt_var_.add_dimension("x", 1, 1, 0);
  ts_var_.add_dimension("x", 1, 1, 0);

  dx_var_.has_time_dimension(false);
  dy_var_.has_time_dimension(false);
  dz_var_.has_time_dimension(false);
  dt_var_.has_time_dimension(false);
  ts_var_.has_time_dimension(false);

  dx_var_.set_datatype(GRID_MPI_TYPE);
  dy_var_.set_datatype(GRID_MPI_TYPE);
  dz_var_.set_datatype(GRID_MPI_TYPE);
  dt_var_.set_datatype(GRID_MPI_TYPE);
  ts_var_.set_datatype(MPI_UNSIGNED);

  dx_var_.set_ptr(&dx_);
  dy_var_.set_ptr(&dy_);
  dz_var_.set_ptr(&dz_);
  dt_var_.set_ptr(&dt_);
  ts_var_.set_ptr(&ts_);

  pre_vars_["dx"] = &dx_var_;
  pre_vars_["dy"] = &dy_var_;
  pre_vars_["dz"] = &dz_var_;
  pre_vars_["dt"] = &dt_var_;
}

void PropertiesResult::deinit()
{}

void PropertiesResult::calculate_pre_result(const Grid &grid)
{
  dx_ = grid.get_deltax();
  dy_ = grid.get_deltay();
  dz_ = grid.get_deltaz();
  dt_ = grid.get_deltat();

  ts_ = 0; // no access to this...

  const GridInfo& gi = grid.get_grid_info();
  dimx_ = gi.global_dimx_;
  dimy_ = gi.global_dimy_;
  dimz_ = gi.global_dimz_;
}

ostream& PropertiesResult::to_string(ostream &os) const
{
  return os << "PropertiesResult, grid and time deltas, etc \n";
}
