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

#include "../MetaFDTD.hh"

using namespace boost::python;


void export_fdtd()
{
  class_<FDTD, boost::noncopyable>("FDTD")
    .def("set_grid_size", &FDTD::set_grid_size)
    .def("set_grid_centre", &FDTD::set_grid_centre)
    .def("get_grid_size", &FDTD::get_grid_size)
    .def("get_grid_centre", &FDTD::get_grid_centre)
    .def("set_grid_deltas", &FDTD::set_grid_deltas)
    .def("set_time_delta", &FDTD::set_time_delta)
    .def("set_boundary", &FDTD::set_boundary)
    .def("load_materials", &FDTD::load_materials)
    .def("add_excitation", &FDTD::add_excitation)
    .def("add_result", &FDTD::add_result)
    .def("add_object", &FDTD::add_object)
    .def("add_datawriter", &FDTD::add_datawriter)
    .def("map_result_to_datawriter", &FDTD::map_result_to_datawriter)
    .def("run", &FDTD::run)
    .def("set_time_steps", &FDTD::set_time_steps)
    .def("get_time_steps", &FDTD::get_time_steps)
    .def("get_time_delta", &FDTD::get_time_delta)

    .def("get_x_cells", &FDTD::get_num_x_cells)
    .def("get_y_cells", &FDTD::get_num_y_cells)
    .def("get_z_cells", &FDTD::get_num_z_cells)

    .def("get_grid", &FDTD::get_grid)
    .def("set_decomp_alg", &FDTD::set_decomp_alg)
    .def("get_decomp_alg", &FDTD::get_decomp_alg)
    .add_property("decomp_alg", &FDTD::get_decomp_alg, 
                  &FDTD::set_decomp_alg)
    ;
  
  enum_<MetaType>("MetaType")
    .value("METAFDTD_ONE", METAFDTD_ONE)
    .value("METAFDTD_THREE", METAFDTD_THREE)
    //.value("METAFDTD_SGI_ORIGIN", METAFDTD_SGI_ORIGIN)
    .export_values()
    ;

  class_<MetaFDTD, bases<FDTD>, boost::noncopyable>("MetaFDTD")
    .def("set_metatype", &MetaFDTD::set_metatype)
    .def("get_metatype", &MetaFDTD::get_metatype)
    .add_property("metatype", &MetaFDTD::get_metatype,
                  &MetaFDTD::set_metatype)
    ;

}
