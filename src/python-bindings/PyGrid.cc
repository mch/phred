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

#include <boost/python.hpp>

#include "../Grid.hh"
#include "../FreqGrid.hh"

using namespace boost::python;

// class GridWrap : public Grid
// {
//   PyObject* self;

//   GridWrap(PyObject* self_)
//     : self(self_) {}

//   void apply_boundaries() 
//   { 
//     return call_apply_boundaries<void>(self, "apply_boundaries"); 
//   }

//   void load_materials() 
//   { 
//     return call_load_materials<void>(self, "load_materials"); 
//   }

//   void free_material() 
//   { 
//     return call_free_material<void>(self, "free_material"); 
//   }

//   void setup_grid() 
//   { 
//     return call_setup_grid<void>(self, "setup_grid"); 
//   }

//   void alloc_grid() 
//   { 
//     return call_alloc_grid<void>(self, "alloc_grid"); 
//   }

//   void free_grid() 
//   { 
//     return call_free_grid<void>(self, "free_grid"); 
//   }

//   void update_ex() 
//   { return call_update_ex<void>(self, "update_ex"); }

//   void update_ey() 
//   { return call_update_ey<void>(self, "update_ey"); }

//   void update_ez() 
//   { return call_update_ez<void>(self, "update_ez"); }

//   void update_hx() 
//   { return call_update_hx<void>(self, "update_hx"); }

//   void update_hy() 
//   { return call_update_hy<void>(self, "update_hy"); }

//   void update_hz() 
//   { return call_update_hz<void>(self, "update_hz"); }
// };

void export_grids()
{
  class_<Grid>("Grid");
  class_<FreqGrid, bases<Grid>, boost::noncopyable>("FreqGrid");
}
