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

#include "../Types.hh"

using namespace boost::python;

BOOST_PYTHON_MODULE(types)
{
  enum_<BoundaryCondition>("BoundaryCondition")
    .value("UNKNOWN", UNKNOWN)
    .value("SUBDOMAIN", SUBDOMAIN)
    .value("EWALL", EWALL)
    .value("HWALL", HWALL)
    .value("PML", PML)
    .export_values()
    ;

  enum_<MaterialType>("MaterialType")
    .value("PERF_COND", PERF_COND)
    .value("NON_PERMEABLE", NON_PERMEABLE)
    .value("CONDUCTIVE", CONDUCTIVE)
    .value("DEBYE", DEBYE)
    .value("LORENTZ", LORENTZ)
    .value("DRUDE", DRUDE)
    .export_values()
    ;

  enum_<Face>("Face")
    .value("FRONT", FRONT)
    .value("BACK", BACK)
    .value("LEFT", LEFT)
    .value("RIGHT", RIGHT)
    .value("BOTTOM", BOTTOM)
    .value("TOP", TOP)
    .export_values()
    ;

  enum_<FieldComponent>("FieldComponent")
    .value("EX", FC_EX)
    .value("EY", FC_EY)
    .value("EZ", FC_EZ)
    .value("HX", FC_HX)
    .value("HY", FC_HY)
    .value("HZ", FC_HZ)
    .export_values()
    ;

  enum_<FieldType>("FieldType")
    .value("E", E)
    .value("H", H)
    .value("BOTH", BOTH)
    .export_values()
    ;

  class_<region_t>("region")
    .def_readwrite("xmin", &region_t::xmin)
    .def_readwrite("ymin", &region_t::ymin)
    .def_readwrite("zmin", &region_t::zmin)
    .def_readwrite("xmax", &region_t::xmax)
    .def_readwrite("ymax", &region_t::ymax)
    .def_readwrite("zmax", &region_t::zmax)
    ;

  class_<point_t>("point")
    .def_readwrite("x", &point_t::x)
    .def_readwrite("y", &point_t::y)
    .def_readwrite("z", &point_t::z)
    ;
}
