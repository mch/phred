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

#include "PlaneResult.hh"
#include "../Globals.hh"

PlaneResult::PlaneResult()
  : face_(FRONT), field_(FC_EY), average_(false), 
    avg_data_(0), have_data_(true)
{
  variables_["Plane"] = &var_;
}

PlaneResult::~PlaneResult()
{
  ////deinit(); // deinit() is not safe to call multiple times right now.
}

void PlaneResult::calculate_result(const Grid &grid, 
                                   unsigned int time_step)
{
  if (result_time(time_step) && have_data_)
  {
    var_.set_num(1);

    if (average_)
    {
      int idx = 0;
      switch(face_)
      {
      case FRONT:
        for (int j = (*region_).ymin(); j <= (*region_).ymax(); j++)
        {
          for (int k = (*region_).zmin(); k <= (*region_).zmax(); k++)
          {
            avg_data_ = 0;
            idx++;
          }
        }
        break;
      }
    }

  } else {
    var_.set_num(0);
  }
}

void PlaneResult::init(const Grid &grid)
{
  // The the measurement plane is parallel to all subdomain cut
  // planes, then this PlaneResult can only return data if the point
  // is located in the local subdomain.

  cells_ = grid.get_cellset(*box_);

  region_ = cells_->get_local_block();
  shared_ptr<Block> global_b = cells_->get_global_block();

  have_data_ = (*region_).has_face_data(face_);
  unsigned int sz = 0;
  grid_point gp;

  switch (face_) 
  {
  case FRONT:
  case BACK:
    if (face_ == BACK)
      gp.x = (*region_).xmin();
    else
      gp.x = (*region_).xmax();    
    
    gp.y = (*region_).ymin();
    gp.z = (*region_).zmin();

    // DIMENSION STARTS HAVE TO CHANGE TOO!
    var_.add_dimension("y", (*region_).ylen(), (*global_b).ylen(), 
                       (*region_).yoffset());
    var_.add_dimension("z", (*region_).zlen(), (*global_b).zlen(), 
                       (*region_).zoffset());

    sz = (*region_).ylen() * (*region_).zlen();


    if (!average_)
    {
      MPI_Type_hvector((*region_).ylen(), (*region_).zlen(), 
                       sizeof(field_t) * grid.get_ldz_sd(), 
                       GRID_MPI_TYPE, &datatype_);
    }
    break;

  case TOP:
  case BOTTOM:
    gp.x = (*region_).xmin();
    gp.y = (*region_).ymin();
    
    if (face_ == BOTTOM)
      gp.z = (*region_).zmin();
    else
      gp.z = (*region_).zmax();    
    
    var_.add_dimension("x", (*region_).xlen(), (*global_b).xlen(), 
                       (*region_).xoffset());
    var_.add_dimension("y", (*region_).ylen(), (*global_b).ylen(), 
                       (*region_).yoffset());

    sz = (*region_).ylen() * (*region_).xlen();

    if (!average_)
    {
      MPI_Datatype y_vector;
      MPI_Type_vector((*region_).ylen(), 1, grid.get_ldz_sd(), 
                      GRID_MPI_TYPE, &y_vector);

      MPI_Type_hvector((*region_).xlen(), 1, 
                       sizeof(field_t) * grid.get_ldz_sd() 
                       * grid.get_ldy_sd(), 
                       y_vector, &datatype_);
    }
    break;

  case LEFT:
  case RIGHT:
    gp.x = (*region_).xmin();
    
    if (face_ == LEFT)
      gp.y = (*region_).ymin();
    else
      gp.y = (*region_).ymax();    

    gp.z = (*region_).zmin();

    // 2005-02-12, MCH: Flipped x and z for a little test.... Seems to
    // have fixed the problem.
    var_.add_dimension("x", (*region_).xlen(), (*global_b).xlen(), 
                       (*region_).xoffset());
    var_.add_dimension("z", (*region_).zlen(), (*global_b).zlen(), 
                       (*region_).zoffset());

    sz = (*region_).xlen() * (*region_).zlen();

    if (!average_)
    {
      MPI_Type_vector((*region_).xlen(), (*region_).zlen(), 
                      grid.get_ldy_sd() * grid.get_ldz_sd(), 
                      GRID_MPI_TYPE, &datatype_);
    }

    break;
  }

  if (average_)
  {
    MPI_Type_contiguous(sz, GRID_MPI_TYPE, &datatype_);
  }

  MPI_Type_commit(&datatype_);

  var_.set_name(base_name_);

  if (average_)
  {
    avg_data_ = new field_t[sz];
    memset(avg_data_, 0, sz * sizeof(field_t));
    var_.set_ptr(avg_data_);
  } 
  else
  {
    var_.set_datatype(datatype_);

    var_.set_ptr(grid.get_pointer(gp, field_));
  }
}

void PlaneResult::deinit()
{
  if (datatype_)
    MPI_Type_free(&datatype_);
}

ostream& PlaneResult::to_string(ostream &os) const
{
  //FieldComponent f = field_;
  //Face face = face_;
  os << "PlaneResult returning " << field_component_string(field_)
     << " data on a face parallel to the " 
     << face_string(face_) << ". In grid coordinates, this plane is "
     << *region_ << endl;

  return os;
}
