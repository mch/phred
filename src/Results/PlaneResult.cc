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
#include "../Globals.hh"

PlaneResult::PlaneResult()
  : face_(FRONT), field_(FC_EY), have_data_(true)
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
  if (result_time(time_step) && have_data_)
  {
    var_.set_num(1);
  } else {
    var_.set_num(0);
  }

  return variables_;
}

void PlaneResult::init(const Grid &grid)
{
  // The the measurement plane is parallel to all subdomain cut
  // planes, then this PlaneResult can only return data if the point
  // is located in the local subdomain.

//   cerr << "PlaneResult global point location: " << plane_.x 
//        << " x " << plane_.y << " x " << plane_.z << endl;
//   grid_point l_point = grid.global_to_local(plane_);
//   cerr << "PlaneResult local point location: " << l_point.x 
//        << " x " << l_point.y << " x " << l_point.z << endl;

  bool x_sd = false;
  bool y_sd = false;
  bool z_sd = false;

  if (grid.get_ldx() != grid.get_gdx())
  {
    x_sd = true;
  }

  if (grid.get_ldy() != grid.get_gdy())
  {
    y_sd = true;
  }

  if (grid.get_ldz() != grid.get_gdz())
  {
    z_sd = true;
  }

  //cerr << MPI_RANK << ": ";
 
  switch (face_) 
  {
  case FRONT:
  case BACK:
    // DIMENSION STARTS HAVE TO CHANGE TOO!
    var_.add_dimension("y", grid.get_ldy(), grid.get_gdy(), grid.get_lsy_ol());
    var_.add_dimension("z", grid.get_ldz(), grid.get_gdz(), grid.get_lsz_ol());

    cerr << "X subdomin split? " << (x_sd ? "yup":"nope") << " plane_.x: "
         << plane_.x << ", min: " << grid.get_lsx() 
         << ", max: " << grid.get_lsx() + grid.get_ldx() << endl;
     
    if (x_sd && (plane_.x > grid.get_ldx() - 1 + grid.get_lsx() 
                 || plane_.x < grid.get_lsx()))
      have_data_ = false;

    // ERROR: This is only contiguous if the overlap IS INCLUDED!
    // That's not what we want for results!
    MPI_Type_contiguous(grid.get_ldz() * grid.get_ldy(), 
                        GRID_MPI_TYPE, &datatype_);

    //MPI_Type_vector(grid.get_ldy(), grid.get_ldz(), 
    //                grid.get_ldz_sd() - grid.get_ldz(), 
    //                GRID_MPI_TYPE, &datatype_);

    break;

  case TOP:
  case BOTTOM:
    var_.add_dimension("x", grid.get_ldx(), grid.get_gdx(), grid.get_lsx_ol());
    var_.add_dimension("y", grid.get_ldy(), grid.get_gdy(), grid.get_lsy_ol());

    MPI_Datatype y_vector;
    MPI_Type_vector(grid.get_ldy(), 1, grid.get_ldz_sd(), 
                    GRID_MPI_TYPE, &y_vector);

    MPI_Type_hvector(grid.get_ldx(), 1, 
                     sizeof(field_t) * grid.get_ldz_sd() * grid.get_ldy_sd(), 
                     y_vector, &datatype_);

     cerr << "Z subdomin split? " << (z_sd ? "yup":"nope") << " plane_.z: "
          << plane_.z << ", min: " << grid.get_lsz() 
          << ", max: " << grid.get_lsz() + grid.get_ldz() << endl;

    if (z_sd && (plane_.z > grid.get_ldz() - 1 + grid.get_lsz() 
                 || plane_.z < grid.get_lsz()))
      have_data_ = false;

    break;

  case LEFT:
  case RIGHT:
    var_.add_dimension("x", grid.get_ldx(), grid.get_gdx(), grid.get_lsx_ol());
    var_.add_dimension("z", grid.get_ldz(), grid.get_gdz(), grid.get_lsz_ol());

    cerr << "Y subdomin split? " << (y_sd ? "yup":"nope") << " plane_.y: "
         << plane_.y << ", min: " << grid.get_lsy() 
         << ", max: " << grid.get_lsy() + grid.get_ldy() << endl;

    if (y_sd && (plane_.y > grid.get_ldy() + grid.get_lsy() 
                 || plane_.y < grid.get_lsy()))
      have_data_ = false;

    MPI_Type_vector(grid.get_ldx(), grid.get_ldz(), 
                    grid.get_ldy() * grid.get_ldz(), 
                    GRID_MPI_TYPE, &datatype_);
    break;
  }

  cerr << "Have data? " << (have_data_ ? "yes" : "no") << endl;

  MPI_Type_commit(&datatype_);

  var_.set_name(base_name_);
  var_.set_datatype(datatype_);
  var_.set_ptr(grid.get_face_start(face_, field_, plane_));

}

void PlaneResult::deinit()
{
  MPI_Type_free(&datatype_);
}
