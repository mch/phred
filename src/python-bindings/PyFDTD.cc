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

#include "../FDTD.hh"

using namespace boost::python;


void export_fdtd()
{
  class_<FDTD>("FDTD")
    .def("set_grid_size", &FDTD::set_grid_size)
    .def("set_grid_deltas", &FDTD::set_grid_deltas)
    .def("set_boundary", &FDTD::set_boundary)
    .def("load_materials", &FDTD::load_materials)
    .def("add_excitation", &FDTD::add_excitation)
    .def("add_result", &FDTD::add_result)
    .def("add_geometry", &FDTD::add_geometry)
    .def("add_datawriter", &FDTD::add_datawriter)
    .def("map_result_to_datawriter", &FDTD::map_result_to_datawriter)
    .def("run", &FDTD::run)
    .def("set_time_steps", &FDTD::set_time_steps)
    ;

  
}
