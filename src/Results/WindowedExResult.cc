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

#include "WindowedExResult.hh"
#include "../Globals.hh"

WindowedExResult::WindowedExResult(shared_ptr<WindowedExcitation> wex)
  : wex_(wex)
{}

void WindowedExResult::init(const Grid &grid)
{
  var_.reset();
  pre_vars_["WindowResult"] = &var_;

  shared_ptr<CSGBox> box = wex_->get_region();
  shared_ptr<CellSet> cells = grid.get_cellset(*(box.get()));
  gregion_ = cells->get_global_block();

  var_.add_dimension("x", gregion_->xlen(), gregion_->xlen(), 0);
  var_.add_dimension("y", gregion_->ylen(), gregion_->ylen(), 0);
  var_.add_dimension("z", gregion_->zlen(), gregion_->zlen(), 0);
  var_.has_time_dimension(false);

  int sz = gregion_->xlen() * gregion_->ylen() * gregion_->zlen();

  wnd_ = new float[sz];

  var_.set_ptr(wnd_);
  var_.set_num(sz);
  //var_.set_element_type(MPI_FLOAT);
  var_.set_name(base_name_);
}

void WindowedExResult::deinit()
{
  if (wnd_)
  {
    delete[] wnd_;
    wnd_ = 0;
  }
}

void WindowedExResult::calculate_pre_result(const Grid &grid)
{
  if (MPI_RANK == 0)
  {
    delta_t dx, dy, dz;
    dx = grid.get_deltax();
    dy = grid.get_deltay();
    dz = grid.get_deltaz();
    
    float fx, fy, fz, xmin, ymin, zmin;

    int idx = 0;
    
    shared_ptr<CSGBox> box = wex_->get_region();

    point centre = box->get_centre();
    point size = box->get_size();
  
    xmin = centre.x - size.x / 2;
    ymin = centre.y - size.y / 2;
    zmin = centre.z - size.z / 2;

    fx = xmin;
    for(loop_idx_t i = (*gregion_).xmin(); i <= (*gregion_).xmax(); i++, 
          fx += dx)
    {
      fy = ymin;
      for (loop_idx_t j = (*gregion_).ymin(); j <= (*gregion_).ymax(); j++,
             fy += dy)
      {
        fz = zmin;
        for (loop_idx_t k = (*gregion_).zmin(); k <= (*gregion_).zmax(); k++,
               fz += dz)
        {
          wnd_[idx++] = wex_->window(fx, fy, fz);
        }
      }
    }

  } else {
    var_.set_num(0);
  }
}


ostream& WindowedExResult::to_string(ostream &os) const
{
  return os << "WindowedExResult returning an evaluated windowing function..."
            << endl;
}
