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

#include "BlockResult.hh"
#include "../Globals.hh"
#include "../Contiguous.hh"

#include <math.h>

BlockResult::BlockResult()
  : field_comp_(FC_EY), init_(false), field_data_(0)
{
  variables_["block"] = &var_;
}

BlockResult::BlockResult(region_t r, FieldComponent field_comp)
  : region_(r), field_comp_(field_comp), init_(false), field_data_(0)
{
  variables_["block"] = &var_;
}

BlockResult::~BlockResult()
{
  deinit();
}

void BlockResult::deinit()
{
  if (field_data_) {
    delete[] field_data_;
    field_data_ = 0;
  }
}

// void BlockResult::init(const Grid &grid)
// {
//   MPI_Datatype temp;
//   int sizes[3];
//   int subsizes[3];
//   int g_subsizes[3];
//   int starts[3];
  
//   g_subsizes[0] = region_.xmax - region_.xmin;
//   g_subsizes[1] = region_.ymax - region_.ymin;
//   g_subsizes[2] = region_.zmax - region_.zmin;

//   // Setup (convert to local)
//   region_ = grid.global_to_local(region_);
//   sizes[0] = grid.get_ldx();
//   sizes[1] = grid.get_ldy();
//   sizes[2] = grid.get_ldz();

//   subsizes[0] = region_.xmax - region_.xmin;
//   subsizes[1] = region_.ymax - region_.ymin;
//   subsizes[2] = region_.zmax - region_.zmin;

//   starts[0] = region_.xmin;
//   starts[1] = region_.ymin;
//   starts[2] = region_.zmin;

//   // Create
//   MPI_Type_create_subarray(3, sizes, subsizes, starts, 1, 
//                            GRID_MPI_TYPE, &temp);
//   MPI_Type_commit(&temp);

//   var_.set_name(base_name_);
//   var_.set_datatype(temp);
//   var_.set_num(0);
//   var_.set_ptr(const_cast<field_t *>(grid.get_pointer(grid_point(region_.xmin, 
//                                                               region_.ymin, 
//                                                               region_.zmin), 
//                                                       field_comp_)));
//   var_.add_dimension("x", subsizes[0], g_subsizes[0], starts[0]);
//   var_.add_dimension("y", subsizes[1], g_subsizes[1], starts[1]);
//   var_.add_dimension("z", subsizes[2], g_subsizes[2], starts[2]);
  
//   init_ = true; 
// }

void BlockResult::init(const Grid &grid)
{
//   if (MPI_SIZE != 1)
//     throw DataWriterException("BlockResult is only available when running on a single node for now.");
  
  if (MPI_SIZE == 1)
  {
    MPI_Type_contiguous(grid.get_ldx() * grid.get_ldy() * grid.get_ldz(), 
                        GRID_MPI_TYPE, &datatype_);
    MPI_Type_commit(&datatype_);
  } 
  else
  {
    unsigned int *dsize = new unsigned int[3];
    unsigned int *coord = new unsigned int[3];
    unsigned int *lens = new unsigned int[3];
    unsigned int ndisps = 0;
    int *displacements = 0;
    const GridInfo &info = grid.get_grid_info();

    dsize[0] = info.dimx_;
    dsize[1] = info.dimy_;
    dsize[2] = info.dimz_;

    coord[0] = info.start_x_no_sd_ - info.start_x_;
    coord[1] = info.start_y_no_sd_ - info.start_y_;
    coord[2] = info.start_z_no_sd_ - info.start_z_;

    lens[0] = info.dimx_no_sd_;
    lens[1] = info.dimy_no_sd_;
    lens[2] = info.dimz_no_sd_;

    Contiguous::compute_displacements(3, dsize, coord, lens, &ndisps, 
                                      &displacements);

    if (ndisps > 0)
    {
      int *contig_lens = new int[ndisps];

      for (int i = 0; i < ndisps; i++)
        contig_lens[i] = static_cast<int>(info.dimz_no_sd_);

      MPI_Type_indexed(ndisps, contig_lens, displacements, 
                       GRID_MPI_TYPE, &datatype_);
      MPI_Type_commit(&datatype_);
    } else {
      cerr << "BlockResult error: number of displacements is zero. ";
    }
  }

  var_.add_dimension("x", grid.get_ldx(), grid.get_gdx(), grid.get_lsx());
  var_.add_dimension("y", grid.get_ldy(), grid.get_gdy(), grid.get_lsy());
  var_.add_dimension("z", grid.get_ldz(), grid.get_gdz(), grid.get_lsz());
  
  var_.set_name(base_name_);
  var_.set_datatype(datatype_);
  
  if (field_comp_ == FC_E || field_comp_ == FC_H)
  {
    field_data_ = new field_t[grid.get_ldx() 
                              * grid.get_ldy() 
                              * grid.get_ldz()];
    var_.set_ptr(field_data_);
  } 
  else
  {
    var_.set_ptr(const_cast<field_t *>
                 (grid.get_pointer(grid_point(0,0,0), field_comp_)));
  }

  init_ = true;
}

map<string, Variable *> &BlockResult::get_result(const Grid &grid, 
                                         unsigned int time_step)
{
  if (init_ && result_time(time_step))
  {
    if (field_comp_ == FC_E)
    {
      field_t *f = field_data_;
      const field_t *ex = grid.get_pointer(grid_point(0,0,0), FC_EX);
      const field_t *ey = grid.get_pointer(grid_point(0,0,0), FC_EY);
      const field_t *ez = grid.get_pointer(grid_point(0,0,0), FC_EZ);

      for (unsigned int i = 0; i < grid.get_ldx(); i++)
      {
        for (unsigned int j = 0; j < grid.get_ldy(); j++)
        {
          ex = grid.get_pointer(grid_point(i, j, 0), FC_EX);
          ey = grid.get_pointer(grid_point(i, j, 0), FC_EY);
          ez = grid.get_pointer(grid_point(i, j, 0), FC_EZ);
          for (unsigned int k = 0; k < grid.get_ldz(); 
               k++, f++, ex++, ey++, ez++)
          {
            *f = sqrt(pow(*ex, 2) + pow(*ey, 2) + pow(*ez, 2));
//             field_data_[grid.pi(i,j,k)] 
//               = sqrt(pow(grid.get_ex(i,j,k), 2) 
//                      + pow(grid.get_ey(i,j,k), 2) 
//                      + pow(grid.get_ez(i,j,k), 2) );
          }
        }
      }
    }
    else if (field_comp_ == FC_H)
    {
      field_t *f = field_data_;
      const field_t *hx = grid.get_pointer(grid_point(0,0,0), FC_HX);
      const field_t *hy = grid.get_pointer(grid_point(0,0,0), FC_HY);
      const field_t *hz = grid.get_pointer(grid_point(0,0,0), FC_HZ);

      for (unsigned int i = 0; i < grid.get_ldx(); i++)
        for (unsigned int j = 0; j < grid.get_ldy(); j++)
          for (unsigned int k = 0; k < grid.get_ldz(); 
               k++, f++, hx++, hy++, hz++)
            *f = sqrt(pow(*hx, 2) + pow(*hy, 2) + pow(*hz, 2));

//             field_data_[grid.pi(i,j,k)] 
//               = sqrt(pow(grid.get_hx(i,j,k), 2) 
//                      + pow(grid.get_hy(i,j,k), 2) 
//                      + pow(grid.get_hz(i,j,k), 2) );
    }

    var_.set_num(1);
  } 
  else 
  {
    var_.set_num(0);
  }

  return variables_;
}
