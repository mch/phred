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
#include <boost/python/args.hpp>

#include "../MetaFDTD.hh"

using namespace boost::python;


void export_fdtd()
{
  class_<FDTD, boost::noncopyable>("FDTD", "This object holds all "
                                   "information about a simulation.")
    .def("set_grid_size", &FDTD::set_grid_size, 
         "Set the size of the FDTD grid in meters. set_grid_size(x, y, z)")
    .def("set_grid_centre", &FDTD::set_grid_centre,
         "Set the centre of the grid in meters. set_grid_centre(x, y, z)")
    .def("get_grid_size", &FDTD::get_grid_size,
         "Returns the grid size as a point object")
    .def("get_grid_centre", &FDTD::get_grid_centre, 
         "Returns the centre of the grid as a point object")
    .def("set_grid_deltas", &FDTD::set_grid_deltas,
         "Set the grid cell sizes. set_grid_delta(dx, dy, dz)")
    .def("set_grid_material", &FDTD::set_grid_material,
         "Set the default material for the interior of the grid.")
    .def("set_time_delta", &FDTD::set_time_delta,
         "Set the time step size. The time step size is normally "
         "computed from the cell size. set_time_delta(dt)")
    .def("set_boundary", &FDTD::set_boundary,
         "Set a boundary conditon on a face. Face can be one of FRONT, "
         "BACK, LEFT, RIGHT, TOP, or BOTTOM, and the obj can be a boundary "
         "condition object such as UPml, Ewall, Hwall, or Periodic. "
         "set_boundary(Face, obj)")
    .def("load_materials", &FDTD::load_materials,
         "Load a material library object. load_materials(matlib)")
    .def("add_excitation", &FDTD::add_excitation,
         "Add an excitation object to the simulation. "
         "add_excitation(name, obj)")
    .def("add_result", 
         (void(FDTD::*)(const char *name, shared_ptr<Result>))
         &FDTD::add_result, 
         "Add a result object. add_result(name, obj)")
    .def("add_result_and_map", 
         (void(FDTD::*)(const char *name, shared_ptr<Result>,
                        const char *dw))&FDTD::add_result, 
         "Add a result object and map it directly to a datawriter. "
         "add_result_and_map(name, obj, dw_name)")
    .def("add_object", &FDTD::add_object, 
         "Add a geometry object to the simulation. "
         "add_object(material_name, obj)")
    .def("add_datawriter", &FDTD::add_datawriter,
         "Add a datawriter. add_datawriter(name, obj)")
    .def("map_result_to_datawriter", &FDTD::map_result_to_datawriter, 
         "Map a result to a datawriter. "
         "map_result_to_datawriter(result_name, datawriter_name)")
    .def("run", &FDTD::run, "Run the simulation.")
    .def("set_time_steps", &FDTD::set_time_steps, 
         "Set the number of time steps to run the simulation for.")
    .def("get_time_steps", &FDTD::get_time_steps,
         "Returns the number of time steps the simulation will run for.")
    .def("get_time_delta", &FDTD::get_time_delta,
         "Returns the time step size.")

    .def("get_x_cells", &FDTD::get_num_x_cells, 
         "Returns the number of cells in the domain along the x axis.")
    .def("get_y_cells", &FDTD::get_num_y_cells, 
         "Returns the number of cells in the domain along the y axis.")
    .def("get_z_cells", &FDTD::get_num_z_cells, 
         "Returns the number of cells in the domain along the z axis.")

    .def("get_grid", &FDTD::get_grid, 
         "Returns the grid object. You should never have to use this.")
    .add_property("dt_scale", &FDTD::get_dt_scale, 
                  &FDTD::set_dt_scale /*, 
                  "A scale factor used in the calculation of the time "
                  "step size. Range is 0 < sf < 1. Defaults to 0.9." */);

  // DEPRECATED:
//     .def("set_decomp_alg", &FDTD::set_decomp_alg)
//     .def("get_decomp_alg", &FDTD::get_decomp_alg)
//     .add_property("decomp_alg", &FDTD::get_decomp_alg, 
//                   &FDTD::set_decomp_alg)

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
