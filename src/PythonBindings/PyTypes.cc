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

#include "../Types.hh"

using namespace boost::python;

void export_types()
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
    .value("DIELECTRIC", DIELECTRIC)
    .value("LOSSY", LOSSY)
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
    .value("FC_EX", FC_EX)
    .value("FC_EY", FC_EY)
    .value("FC_EZ", FC_EZ)
    .value("FC_HX", FC_HX)
    .value("FC_HY", FC_HY)
    .value("FC_HZ", FC_HZ)
    .value("FC_E", FC_E)
    .value("FC_H", FC_H)
    .export_values()
    ;

  enum_<FieldType>("FieldType")
    .value("E_FIELD", E)
    .value("H_FIELD", H)
    .value("BOTH", BOTH)
    .export_values()
    ;

  enum_<DomainDecompAlg>("DomainDecompAlg")
    .value("DDA_UNDEFINED", DDA_UNDEFINED)
    .value("DDA_SIMPLE", DDA_SIMPLE)
    .value("DDA_MPICART", DDA_MPICART)
    .value("DDA_STRIPING", DDA_STRIPING)
    .export_values();

  class_<region_t>("region")
    .def_readwrite("xmin", &region_t::xmin)
    .def_readwrite("ymin", &region_t::ymin)
    .def_readwrite("zmin", &region_t::zmin)
    .def_readwrite("xmax", &region_t::xmax)
    .def_readwrite("ymax", &region_t::ymax)
    .def_readwrite("zmax", &region_t::zmax)
    ;

  class_<grid_point>("grid_point")
    .def(init<unsigned int, unsigned int, unsigned int>())
    .def_readwrite("x", &grid_point::x)
    .def_readwrite("y", &grid_point::y)
    .def_readwrite("z", &grid_point::z)
    ;

  class_<point>("point")
    .def(init<float, float, float>())
    .def_readwrite("x", &point::x)
    .def_readwrite("y", &point::y)
    .def_readwrite("z", &point::z)
    ;

}
