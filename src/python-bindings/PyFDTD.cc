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

/**
 * This wrapper handles reference counting and stuff for Python
   objects that are passed into it for storgage and use in the
   simulation.
 */
class FDTDWrap : public FDTD
{
private:
  PyObject* self_;

  // These objects are here simply to increment the reference count to
  // the Python objects so that they don't go out of scope.
  map<string, object> py_excitations_;
  map<string, object> py_results_;
  map<string, object> py_datawriters_;
  //map<string, object> py_geometry_;
  vector<object> py_geometry_;
  object py_mlib_;
  object py_bound_[6];

public:
  FDTDWrap(PyObject* self)
    : self_(self) {}

  void set_boundary(Face face, object bc)
  {
    py_bound_[face] = bc;
    FDTD::set_boundary(face, extract<BoundaryCond*>(bc));
  }

  void load_materials(object matlib)
  {
    py_mlib_ = matlib;
    FDTD::load_materials(extract<MaterialLib&>(matlib));
  }

  void add_excitation(const char *name, object ex)
  {
    py_excitations_[name] = ex;
    FDTD::add_excitation(name, extract<Excitation*>(ex));
  }

  void add_result(const char *name, object r)
  {
    py_results_[name] = r;
    FDTD::add_result(name, extract<Result*>(r));
  }

  void add_datawriter(const char *name, object dw)
  {
    py_datawriters_[name] = dw;
    FDTD::add_datawriter(name, extract<DataWriter*>(dw));
  }

  //void add_geometry(const char *name, object g);
  void add_geometry(object g)
  {
    py_geometry_.push_back(g);
    FDTD::add_geometry(extract<Geometry*>(g));
  }
};

void export_fdtd()
{
  class_<FDTD, FDTDWrap, boost::noncopyable>("FDTD")
    .def("set_grid_size", &FDTD::set_grid_size)
    .def("set_grid_deltas", &FDTD::set_grid_deltas)
    .def("set_boundary", &FDTDWrap::set_boundary)
    .def("load_materials", &FDTDWrap::load_materials)
    .def("add_excitation", &FDTDWrap::add_excitation)
    .def("add_result", &FDTDWrap::add_result)
    .def("add_geometry", &FDTDWrap::add_geometry)
    .def("add_datawriter", &FDTDWrap::add_datawriter)
    .def("map_result_to_datawriter", &FDTD::map_result_to_datawriter)
    .def("run", &FDTD::run)
    .def("set_time_steps", &FDTD::set_time_steps)
    .def("get_time_delta", &FDTD::get_time_delta)
    ;

  
}
