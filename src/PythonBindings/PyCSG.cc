/* 
   Phred - Phred is a parallel finite difference time domain
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

using namespace boost::python;

#include "../CSG/CSGBox.hh"
#include "../CSG/CSGSphere.hh"

void export_csg() 
{
  class_<CSGObject>("CSGObject", "Base class for all CSGObjects.", 
                    no_init)
    .def("is_point_inside", &CSGObject::is_point_inside)
    ;
  
  class_<CSGPrimitive, bases<CSGObject> >("CSGPrimitive",
                       "Base class for CSG primitives.",
                       no_init)
    .def("is_point_inside", &CSGPrimitive::is_point_inside)
    .def("set_centre", &CSGPrimitive::set_centre)
    .def("get_centre", &CSGPrimitive::get_centre)
    ;

  class_<CSGBox, bases<CSGPrimitive> >("CSGBox", 
                                       "A CSG object representing a box.")
    .def("set_size", &CSGBox::set_size)
    .def("get_size", &CSGBox::get_size)
    ;

  class_<CSGSphere, bases<CSGPrimitive> >("CSGSphere",
                                          "A CSG object representing a sphere.")
    .def("set_radius", &CSGSphere::set_radius)
    .def("get_radius", &CSGSphere::get_radius)
    ;
}

