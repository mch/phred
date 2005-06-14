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

#include <boost/python.hpp>

#include "../Grid.hh"
#include "../FreqGrid.hh"

using namespace boost::python;

void export_grids()
{
  class_<Grid, boost::noncopyable>("Grid", no_init)
    .def("get_define_mode", &Grid::get_define_mode)
    .def("set_define_mode", &Grid::set_define_mode)

    .def("get_deltax", &Grid::get_deltax)
    .def("get_deltay", &Grid::get_deltay)
    .def("get_deltaz", &Grid::get_deltaz)
    .def("get_deltat", &Grid::get_deltat)

    .def("set_ex", &Grid::set_ex)
    .def("set_ey", &Grid::set_ey)
    .def("set_ez", &Grid::set_ez)
    .def("set_hx", &Grid::set_hx)
    .def("set_hy", &Grid::set_hy)
    .def("set_hz", &Grid::set_hz)

    .def("get_ex", &Grid::get_ex)
    .def("get_ey", &Grid::get_ey)
    .def("get_ez", &Grid::get_ez)
    .def("get_hx", &Grid::get_hx)
    .def("get_hy", &Grid::get_hy)
    .def("get_hz", &Grid::get_hz)
    
    .def("set_material", &Grid::set_material)
    .def("get_material", &Grid::get_material)

    .def("setup_gird", &Grid::setup_grid)
    .def("load_materials", &Grid::load_materials)
    .def("load_geometry", &Grid::load_geometry)
    ;
  class_<FreqGrid, bases<Grid>, boost::noncopyable >("FreqGrid", no_init);
}
