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
#include <boost/ref.hpp>

#include "../Excitation.hh"
#include "../WindowedExcitation.hh"
#include "../Gaussm.hh"

using namespace boost::python;

/**
 * Helper function for Python classes derived from Excitation
 */
void call_excite(Excitation& ex, Grid &grid, 
                 unsigned int time_step, FieldType type) 
{ return ex.excite(grid, time_step, type); }

/**
 * Helper function for Python classes derived from SourceFunction
 */
field_t call_source_function(SourceFunction& sf, Grid &grid, 
                             unsigned int time_step) 
{ return sf.source_function(grid, time_step); }

/**
 * This wrapper allows for derived classes built in Python.
 */
class ExcitationWrap : public Excitation
{
  PyObject* self_;

public:
  ExcitationWrap(PyObject* self, SourceFunction *sf)
    :  Excitation(sf), self_(self) {}

  ExcitationWrap(PyObject *self, const Excitation &e)
    : Excitation(e), self_(self)
  {}

  void excite(Grid &grid, unsigned int time_step,
              FieldType type) 
  { call_method<void>(self_, "excite", boost::ref(grid), time_step, type); }

  void default_excite(Grid &grid, unsigned int time_step,
                      FieldType type) 
  { Excitation::excite(grid, time_step, type); }
};

/**
 * Free helper function for calling window() on WindowedExcitation's. 
 */
field_t call_window(WindowedExcitation &w, 
                    region_t r, unsigned int x, 
                    unsigned int y, unsigned int z)
{ w.window(r, x, y, z); }

/**
 * This wrapper allows subclasses of WindowedExcitation written in Python
 */
class WindowedExcitationWrap : public WindowedExcitation
{
private:
  PyObject *self_;

public:
  WindowedExcitationWrap(PyObject *self, SourceFunction *sf)
    : WindowedExcitation(sf)
  {}

  void excite(Grid &grid, unsigned int time_step,
              FieldType type) 
  { return call_method<void>(self_, "excite", time_step, type); }

  void default_excite(Grid &grid, unsigned int time_step,
                      FieldType type) 
  { WindowedExcitation::excite(grid, time_step, type); }  

  field_t window(region_t r, unsigned int x, unsigned int y, 
                 unsigned int z)
  { return call_method<field_t>(self_, "window", x, y, z); }
};

/**
 * This wrapper allows for subclasses of SourceFunction written in
 * Python
 */
class SourceFunctionWrap : public SourceFunction
{
private:
  PyObject *self_;

public:
  SourceFunctionWrap(PyObject *self)
    : self_(self)
  {}

  field_t source_function(const Grid &grid, unsigned int time_step)
  { return call_method<field_t>(self_, "source_function"); }
};

void export_excitations()
{
  //def("call_excite", call_excite);
  //def("call_source_function", call_source_function);
    
  class_<Excitation, ExcitationWrap>("Excitation", "Excitations applied to the FDTD grid", init<SourceFunction *>())
    .def("excite", &ExcitationWrap::excite)
    .def("excite", &ExcitationWrap::default_excite)
    .def("set_polarization", &Excitation::set_polarization)
    .def("set_type", &Excitation::set_type)
    .def("set_soft", &Excitation::set_soft)
    .def("get_soft", &Excitation::get_soft)
    .def("set_region", (void(Excitation::*)(region_t))&Excitation::set_region)
    .def("set_region", (void(Excitation::*)(unsigned int, unsigned int, 
                                            unsigned int, unsigned int,
                                            unsigned int, unsigned int))&Excitation::set_region)
    ;

  class_<WindowedExcitation, WindowedExcitationWrap, bases<Excitation>,
    boost::noncopyable>("WindowedExcitation", 
                        "Excitations that apply a windowing function to the excitation in the FDTD grid", 
                        init<SourceFunction *>())
    .def("excite", &WindowedExcitation::excite, 
         &WindowedExcitationWrap::default_excite)
    ;

  class_<SourceFunction, SourceFunctionWrap, boost::noncopyable>("SourceFunction", "Make derived classes from this to create source functions for excitations")
    .def("source_function", &SourceFunction::source_function)
    ;
  //.def("call_sf", call_sf)

  class_<Gaussm, bases<SourceFunction> >("Gaussm", "Gaussian modulated sine function")
    .def("set_parameters", &Gaussm::set_parameters)
    .def("get_alpha", &Gaussm::get_alpha)
    .def("get_deltaf", &Gaussm::get_deltaf)
    .def("get_f0", &Gaussm::get_f0)
    .def("source_function", &Gaussm::source_function)
    ;
}
