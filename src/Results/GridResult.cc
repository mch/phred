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

#include "GridResult.hh"
#include "../Globals.hh"
#include "../Contiguous.hh"

GridResult::GridResult()
{}

GridResult::~GridResult()
{}

void GridResult::init(const Grid &grid)
{
  material_ids_.has_time_dimension(false);
  pre_vars_["material_ids"] = &material_ids_;
  //variables_["delta_xs"] = &deltaxs_;
  //variables_["delta_ys"] = &deltays_;
  //variables_["delta_zs"] = &deltazs_;

#ifdef DEBUG
  cout << "GridResult is returning a region that is "
       << grid.get_ldx() << " by " << grid.get_ldy()
       << " by " << grid.get_ldz() << ".\nEach element is " 
       << sizeof(mat_idx_t) << " bytes. Total bytes: "
       << sizeof(mat_idx_t) * grid.get_ldx() * grid.get_ldy() 
    * grid.get_ldz() << endl;
#endif

  if (MPI_SIZE == 1)
  {
    MPI_Type_contiguous(grid.get_ldx() * grid.get_ldy() * grid.get_ldz(), 
                        MAT_IDX_MPI_TYPE, &datatype_);
    MPI_Type_commit(&datatype_);
  } 
  else
  {
    int *dsize = new int[3];
    int *coord = new int[3];
    int *lens = new int[3];
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

//     Contiguous::compute_displacements(3, dsize, coord, lens, &ndisps, 
//                                       &displacements);

//     if (ndisps > 0)
//     {
//       int *contig_lens = new int[ndisps];

//       for (int i = 0; i < ndisps; i++)
//         contig_lens[i] = static_cast<int>(info.dimz_no_sd_);

//       MPI_Type_indexed(ndisps, contig_lens, displacements, 
//                        MAT_IDX_MPI_TYPE, &datatype_);
//       MPI_Type_commit(&datatype_);
//     } else {
//       cerr << "BlockResult error: number of displacements is zero. ";
//     }

//     delete[] displacements;

    MPI_Type_create_subarray(3, dsize, lens, coord, MPI_ORDER_C,
                             MAT_IDX_MPI_TYPE, &datatype_);
    MPI_Type_commit(&datatype_);
  }

  material_ids_.set_element_type(MAT_IDX_MPI_TYPE);

  material_ids_.add_dimension("x", grid.get_ldx(), grid.get_gdx(), 
                              grid.get_lsx());
  material_ids_.add_dimension("y", grid.get_ldy(), grid.get_gdy(), 
                              grid.get_lsy());
  material_ids_.add_dimension("z", grid.get_ldz(), grid.get_gdz(), 
                              grid.get_lsz());
  
  material_ids_.set_name(base_name_ + "_material_ids");
  material_ids_.set_datatype(datatype_);
  material_ids_.set_ptr(const_cast<mat_idx_t *>(grid.get_material_ptr(grid_point(0,0,0))));
  
  //deltaxs_.add_dimension

  material_ids_.set_num(1);
  deltaxs_.set_num(0);
  deltays_.set_num(0);
  deltazs_.set_num(0);
}

void GridResult::deinit()
{}

ostream& GridResult::to_string(ostream &os) const
{
  return os << "GridResult; exports material information in the grid.";
}
