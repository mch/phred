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

#include "PlaneResult.hh"

PlaneResult::PlaneResult()
  : face_(FRONT), field_(FC_EY)
{
  plane_.x = 0;
  plane_.y = 0;
  plane_.z = 0;

  variables_["Plane"] = &var_;
}

PlaneResult::~PlaneResult()
{}

map<string, Variable *> &PlaneResult::get_result(const Grid &grid, 
                                                 unsigned int time_step)
{
  if (result_time(time_step))
  {
    var_.set_num(1);
  } else {
    var_.set_num(0);
  }

  return variables_;
}

void PlaneResult::init(const Grid &grid)
{
  switch (face_) 
  {
  case FRONT:
  case BACK:
    // DIMENSION STARTS HAVE TO CHANGE TOO!
    var_.add_dimension("y", grid.get_ldy(), grid.get_gdy(), grid.get_lsy());
    var_.add_dimension("z", grid.get_ldz(), grid.get_gdz(), grid.get_lsz());

    // ERROR: This is only contiguous if the overlap IS INCLUDED!
    // That's not what we want for results!
    MPI_Type_contiguous(grid.get_ldz() * grid.get_ldy(), 
                        GRID_MPI_TYPE, &datatype_);
    break;

  case TOP:
  case BOTTOM:
    var_.add_dimension("x", grid.get_ldx(), grid.get_gdx(), grid.get_lsx());
    var_.add_dimension("y", grid.get_ldy(), grid.get_gdy(), grid.get_lsy());

    MPI_Datatype y_vector;
    MPI_Type_vector(grid.get_ldy(), 1, grid.get_ldz(), 
                    GRID_MPI_TYPE, &y_vector);

    MPI_Type_hvector(grid.get_ldx(), 1, 
                     sizeof(field_t) * grid.get_ldz() * grid.get_ldy(), 
                     y_vector, &datatype_);
    break;

  case LEFT:
  case RIGHT:
    var_.add_dimension("x", grid.get_ldx(), grid.get_gdx(), grid.get_lsx());
    var_.add_dimension("z", grid.get_ldz(), grid.get_gdz(), grid.get_lsz());

    MPI_Type_vector(grid.get_ldx(), grid.get_ldz(), 
                    grid.get_ldy() * grid.get_ldz(), 
                    GRID_MPI_TYPE, &datatype_);
    break;
  }

  MPI_Type_commit(&datatype_);

  var_.set_name(base_name_);
  var_.set_datatype(datatype_);
  var_.set_ptr(grid.get_face_start(face_, field_, plane_));

}

void PlaneResult::deinit()
{
  MPI_Type_free(&datatype_);
}
