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
#include "../CSG/CSGCylinder.hh"
#include "../CSG/CSGDifference.hh"
#include "../CSG/CSGUnion.hh"
#include "../CSG/CSGIntersection.hh"
#include "../CSG/CSGArray.hh"

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

  class_<CSGOperator, bases<CSGObject> >("CSGOperator", no_init)
    ;

  class_<CSGDifference, bases<CSGOperator> >("CSGDifference", 
                                             "Represents a difference between two CSG objects.",
                                             init<shared_ptr<CSGObject>, 
                                             shared_ptr<CSGObject> >())
    ;

  class_<CSGUnion, bases<CSGOperator> >("CSGUnion", 
                                        "Represents a union between two CSG objects.",
                                        init<shared_ptr<CSGObject>, 
                                        shared_ptr<CSGObject> >())
    ;

  class_<CSGIntersection, bases<CSGOperator> >("CSGIntersection", 
                                               "Represents an intersection between two CSG objects.",
                                               init<shared_ptr<CSGObject>, 
                                               shared_ptr<CSGObject> >())
    ;

  class_<CSGArray, bases<CSGObject> >("CSGArray", 
                                      "Represents an array of items",
                                      init<shared_ptr<CSGObject> >())
    .def("set_lengths", &CSGArray::set_lengths)
    .def("set_spacing", &CSGArray::set_spacing)
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

  class_<CSGCylinder, bases<CSGPrimitive> >("CSGCylinder", 
                                            "A CSG cylinder oriented along the z axis")
    .def("set_radius", &CSGCylinder::set_radius)
    .def("get_radius", &CSGCylinder::get_radius)
    .def("set_height", &CSGCylinder::set_height)
    .def("get_height", &CSGCylinder::get_height)
    ;

}

