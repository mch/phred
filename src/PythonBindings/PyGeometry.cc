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
#include "../Box.hh"
#include "../SphereGeom.hh"

using namespace boost::python;

void call_init(Geometry &geom, const Grid &grid)
{ return geom.init(grid); }

void call_deinit(Geometry &geom)
{ return geom.deinit(); }

void call_set_material(Geometry &geom, Grid &grid)
{ return geom.set_material(grid); }

bool call_local_point_inside(Geometry &geom, unsigned int x,
                             unsigned int y, unsigned int z)
{ return geom.local_point_inside(x, y, z); }

/**
 * This is a wrapper which allows Geometry objects to be constructed
 * in Python.
 */
class GeometryWrap : public Geometry
{
private:
  PyObject *self_;

public:
  GeometryWrap(PyObject *self)
    : self_(self) {}

  void init(const Grid &grid)
  { return call_method<void>(self_, "init"); }

  void deinit(const Grid &grid)
  { return call_method<void>(self_, "deinit"); }

  void set_material(Grid &grid)
  { return call_method<void>(self_, "set_material"); }

  bool local_point_inside(unsigned int x, unsigned int y, unsigned int z)
  { return call_method<bool>(self_, "local_point_inside"); } 
};

void export_geometry()
{
  class_<Geometry, GeometryWrap, boost::noncopyable>("Geometry", "Base class for FDTD geometries")
    .add_property("material_id", &Geometry::get_material_id, 
                  &Geometry::set_material_id)
    .def("init", &GeometryWrap::init)
    .def("deinit", &GeometryWrap::deinit)
    .def("set_material", &GeometryWrap::set_material)
    .def("local_point_inside", &GeometryWrap::local_point_inside)
    .def("get_bounding_box", &Geometry::get_bounding_box)
    .def("get_local_bounding_box", &Geometry::get_local_bounding_box)
    ;

  class_<Box, bases<Geometry> >("Box", "A simple box")
    .def(init<region_t>())
    .add_property("material_id", &Box::get_material_id, 
                  &Box::set_material_id)
    .def("init", &Box::init)
    .def("set_material", &Box::set_material)
    .def("local_point_inside", &Box::local_point_inside)
    .def("set_region", (void(Box::*)(unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int))&Box::set_region)
    ;

  class_<Sphere, bases<Geometry> >("Sphere", "A sphere")
    .def(init<point_t, unsigned int>())
    .add_property("radius", &Sphere::get_radius, 
                  &Sphere::set_radius)
    .add_property("centre", &Sphere::get_centre,
                  &Sphere::set_centre)
    .def("init", &Sphere::init)
    .def("set_material", &Sphere::set_material)
    .def("local_point_inside", &Sphere::local_point_inside)
    ;
}
