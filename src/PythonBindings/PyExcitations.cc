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
#include <boost/ref.hpp>

#include "../Excitations/Excitation.hh"
#include "../Excitations/WindowedExcitation.hh"
#include "../Excitations/BartlettExcitation.hh"
#include "../Excitations/WaveguideExcitation.hh"
#include "../Excitations/GaussWindExcitation.hh"
#include "../Signals/Gaussm.hh"
#include "../Signals/ExpSine.hh"

using namespace boost::python;

/**
 * Helper function for Python classes derived from Excitation
 */
void call_excite(Excitation& ex, Grid &grid, 
                 unsigned int time_step, FieldType type) 
{ return ex.excite(grid, time_step, type); }

/**
 * Helper function for Python classes derived from Signal
 */
field_t call_signal_function(Signal& sf, float time) 
{ return sf.signal_function(time); }

/**
 * This wrapper allows for derived classes built in Python.
 */
class ExcitationWrap : public Excitation
{
  PyObject* self_;

public:
  ExcitationWrap(PyObject* self, shared_ptr<Signal> sf)
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
                    float x, float y, float z)
{ return w.window(x, y, z); }

/**
 * This wrapper allows subclasses of WindowedExcitation written in Python
 */
class WindowedExcitationWrap : public WindowedExcitation
{
private:
  PyObject *self_;

public:
  WindowedExcitationWrap(PyObject *self, shared_ptr<Signal> sf)
    : WindowedExcitation(sf)
  {}

  void excite(Grid &grid, unsigned int time_step,
              FieldType type) 
  { return call_method<void>(self_, "excite", boost::ref(grid),
                             time_step, type); }

  void default_excite(Grid &grid, unsigned int time_step,
                      FieldType type) 
  { WindowedExcitation::excite(grid, time_step, type); }  

  field_t window(float x, float y, float z)
  { return call_method<field_t>(self_, "window", x, y, z); }
};

/**
 * This wrapper allows for subclasses of Signal written in
 * Python
 */
class SignalWrap : public Signal
{
private:
  PyObject *self_;

public:
  SignalWrap(PyObject *self)
    : self_(self)
  {}

  field_t signal_function(float time) const
  { return call_method<field_t>(self_, "signal_function", 
                                time); }
};

void export_excitations()
{
  //def("call_excite", call_excite);
  //def("call_signal_function", call_signal_function);
    
  class_<Excitation, ExcitationWrap>("Excitation", "Excitations applied to the FDTD grid", init<shared_ptr<Signal> >())
    .def("excite", &ExcitationWrap::excite)
    .def("excite", &ExcitationWrap::default_excite)
    .def("set_polarization", &Excitation::set_polarization)
    .def("set_type", &Excitation::set_type)
    .def("set_soft", &Excitation::set_soft)
    .def("get_soft", &Excitation::get_soft)
    .def("set_region", &Excitation::set_region)
    ;

  class_<WindowedExcitation, WindowedExcitationWrap, bases<Excitation>,
    boost::noncopyable>("WindowedExcitation", 
                        "Excitations that apply a windowing function to the excitation in the FDTD grid", 
                        init<shared_ptr<Signal> >())
    .def("excite", &WindowedExcitation::excite, 
         &WindowedExcitationWrap::default_excite)
    ;

  class_<Signal, SignalWrap, boost::noncopyable>("Signal", "Make derived classes from this to create signal functions for excitations")
    .def("signal_function", &SignalWrap::signal_function)
    ;
  //.def("call_sf", call_sf)

  class_<Gaussm, bases<Signal> >("Gaussm", "Gaussian modulated sine function")
    .def("set_parameters", &Gaussm::set_parameters)
    .def("get_alpha", &Gaussm::get_alpha)
    .def("get_deltaf", &Gaussm::get_deltaf)
    .def("get_f0", &Gaussm::get_f0)
    .def("length", &Gaussm::length)
    .def("signal_function", &Gaussm::signal_function) // in Signal
    ;

  class_<ExpSine, bases<Signal> >("ExpSine", "Ramping up sine function")
    .def(init<float>())
    .add_property("frequency", &ExpSine::get_frequency, 
                  &ExpSine::set_frequency)
    .add_property("amplitude", &ExpSine::get_amplitude, 
                  &ExpSine::set_amplitude)
    .def("signal_function", &ExpSine::signal_function)
    ;

  class_<BartlettExcitation, bases<WindowedExcitation> >("BartlettExcitation", "Bartlett windowed excitation; an attempt at a plane wave.", init<shared_ptr<Signal> >())
    .def("excite", &BartlettExcitation::excite)
    ;

  class_<WaveguideExcitation, bases<WindowedExcitation> >("WaveguideExcitation", "Waveguide excitation; you know, for those pesky waveguides!", init<shared_ptr<Signal> >())
    .def("excite", &WaveguideExcitation::excite)
    .def("set_mode", &WaveguideExcitation::set_mode)
    ;

  class_<GaussWindExcitation, bases<WindowedExcitation> >("GaussWindow", "Gaussian windowed excitation; approximates a plane wave.", init<shared_ptr<Signal> >())
    .def("excite", &GaussWindExcitation::excite)
    ;
}
