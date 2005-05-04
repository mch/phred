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

#include "BlockResult.hh"
#include "../Globals.hh"
#include "../Contiguous.hh"

#include <cmath>

BlockResult::BlockResult()
  : field_comp_(FC_EY), init_(false), field_data_(0)
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

void BlockResult::init(const Grid &grid)
{
  /* Region must be in out local sub-domain */ 
  shared_ptr<CellSet> cells;
  shared_ptr<Block> global_b;

  if (box_.get())
  {
    cells = grid.get_cellset(*box_);
    region_ = cells->get_local_block();
    global_b = cells->get_global_block();
  } else {
    throw ResultException("BlockResult has no CSGBox defined!");
  }

  var_.add_dimension("x", region_->xlen(), global_b->xlen(), 
                     region_->xoffset());
  var_.add_dimension("y", region_->ylen(), global_b->ylen(), 
                     region_->yoffset());
  var_.add_dimension("z", region_->zlen(), global_b->zlen(), 
                     region_->zoffset());
  
  var_.set_name(base_name_);


  if (field_comp_ == FC_E || field_comp_ == FC_H)
  {
    unsigned int sz = region_->xlen()
      * region_->ylen() 
      * region_->zlen();

    field_data_ = new field_t[sz];
    var_.set_ptr(field_data_);

    MPI_Type_contiguous(sz, GRID_MPI_TYPE, &datatype_);
    MPI_Type_commit(&datatype_);
  } 
  else
  {
    int *dsize = new int[3];
    int *coord = new int[3];
    int *lens = new int[3];
    unsigned int ndisps = 0;
    int *displacements = 0;
    const GridInfo &info = grid.get_grid_info();

    dsize[0] = grid.get_ldx_sd();
    dsize[1] = grid.get_ldy_sd();
    dsize[2] = grid.get_ldz_sd();

    coord[0] = region_->xmin();
    coord[1] = region_->ymin();
    coord[2] = region_->zmin();

    lens[0] = region_->xlen();
    lens[1] = region_->ylen();
    lens[2] = region_->zlen();

//     Contiguous::compute_displacements(3, dsize, coord, lens, &ndisps, 
//                                       &displacements);

//     if (ndisps > 0)
//     {
//       int *contig_lens = new int[ndisps];

//       for (int i = 0; i < ndisps; i++)
//         contig_lens[i] = static_cast<int>(region_->zlen());

//       MPI_Type_indexed(ndisps, contig_lens, displacements, 
//                        GRID_MPI_TYPE, &datatype_);
//       MPI_Type_commit(&datatype_);
//     } else {
//       throw ResultException("BlockResult error: number of displacements is zero. ");
//     }
    MPI_Type_create_subarray(3, dsize, lens, coord, MPI_ORDER_C,
                             GRID_MPI_TYPE, &datatype_);
    MPI_Type_commit(&datatype_);    

    //delete[] displacements;
  
    var_.set_ptr(const_cast<field_t *>
                 (grid.get_pointer(grid_point(0,0,0), 
                                   field_comp_)));
  }

  var_.set_datatype(datatype_);
  
  init_ = true;
}

void BlockResult::calculate_result(const Grid &grid, 
                                   unsigned int time_step)
{
  if (init_ && result_time(time_step))
  {
    if (field_comp_ == FC_E)
    {
      field_t *f = field_data_;
      grid_point gp = grid_point(region_->xmin(),
                                 region_->ymin(),
                                 region_->zmin());

      const field_t *ex = grid.get_pointer(gp, FC_EX);
      const field_t *ey = grid.get_pointer(gp, FC_EY);
      const field_t *ez = grid.get_pointer(gp, FC_EZ);

      for (unsigned int i = region_->xmin(); i <= region_->xmax(); i++)
      {
        for (unsigned int j = region_->ymin(); j <= region_->ymax(); j++)
        {
          ex = grid.get_pointer(grid_point(i, j, region_->zmin()), FC_EX);
          ey = grid.get_pointer(grid_point(i, j, region_->zmin()), FC_EY);
          ez = grid.get_pointer(grid_point(i, j, region_->zmin()), FC_EZ);

          for (unsigned int k = region_->zmin(); k <= region_->zmax(); 
               k++, f++, ex++, ey++, ez++)
          {
            *f = sqrt(*ex * *ex + *ey * *ey + *ez * *ez);
          }
        }
      }
    }
    else if (field_comp_ == FC_H)
    {
      field_t *f = field_data_;
      grid_point gp = grid_point(region_->xmin(),
                                 region_->ymin(),
                                 region_->zmin());

      const field_t *hx = grid.get_pointer(gp, FC_HX);
      const field_t *hy = grid.get_pointer(gp, FC_HY);

      // On the RS/6000 SP, there seems to be a #define for hz!
#undef hz
      const field_t *hz = grid.get_pointer(gp, FC_HZ);

      for (unsigned int i = region_->xmin(); i <= region_->xmax(); i++)
      {
        for (unsigned int j = region_->ymin(); j <= region_->ymax(); j++)
        {
          hx = grid.get_pointer(grid_point(i, j, region_->zmin()), FC_HX);
          hy = grid.get_pointer(grid_point(i, j, region_->zmin()), FC_HY);
          hz = grid.get_pointer(grid_point(i, j, region_->zmin()), FC_HZ);

          for (unsigned int k = region_->zmin(); k <= region_->zmax(); 
               k++, f++, hx++, hy++, hz++)
          {
            *f = sqrt(*hx * *hx + *hy * *hy + *hz * *hz);
          }
        }
      }
    }

    var_.set_num(1);
  } 
  else 
  {
    var_.set_num(0);
  }
}

ostream& BlockResult::to_string(ostream &os) const
{
  os << "BlockResult: returning the ";
  switch(field_comp_)
  {
  case FC_E:
    os << "magnitude of the electric field";
    break;

  case FC_H:
    os << "magnitude of the magnetic field";
    break;

  case FC_EX:
    os << "ex component";
    break;

  case FC_EY:
    os << "ey component";
    break;

  case FC_EZ:
    os << "ez component";
    break;

  case FC_HX:
    os << "hx component";
    break;

  case FC_HY:
    os << "hy component";
    break;

  case FC_HZ:
    os << "hz component";
    break;

  default:
    os << "invalid field";
    break;
  }

  os << " from the following block:" << endl;
  
  os << *region_ << endl;

  return os;
}
