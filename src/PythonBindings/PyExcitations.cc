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
#include <boost/ref.hpp>

#include "../Excitations/Excitation.hh"
#include "../Excitations/WindowedExcitation.hh"
#include "../Excitations/BartlettExcitation.hh"
#include "../Excitations/WaveguideExcitation.hh"
#include "../Excitations/GaussWindExcitation.hh"
#include "../Excitations/PeriodicExcitation.hh"
#include "../Signals/Gaussm.hh"
#include "../Signals/GaussPulse.hh"
#include "../Signals/DiffGaussPulse.hh"
#include "../Signals/ExpSine.hh"

using namespace boost::python;

/**
 * This wrapper allows subclasses of WindowedExcitation written in Python
 */
class WindowedExcitationWrap : public WindowedExcitation //, 
//wrapper<WindowedExcitation>
{
private:
  PyObject *self_;

public:
  WindowedExcitationWrap(PyObject *self, shared_ptr<Signal> sf)
    : self_(self), WindowedExcitation(sf)
  { }

//   WindowedExcitationWrap(shared_ptr<Signal> sf)
//     : WindowedExcitation(sf)
//   { cout << "WindowedExcitationWrap constructor()" << endl; }

  void init(const Grid &grid)
  { 
    WindowedExcitation::init(grid);
    return call_method<void>(self_, "init", boost::ref(grid)); 
    //return this->get_override("init")(grid);

    // if (override init = this->get_override("init"))
//       init(grid);
  }

  field_t window(float x, float y, float z)
  { 
    return call_method<field_t>(self_, "window", x, y, z); 
    //return this->get_override("window")(x, y, z);
  }

  float get_xmin() 
  { return xmin_; }

  float get_ymin() 
  { return ymin_; }

  float get_zmin() 
  { return zmin_; }

  float get_xmax() 
  { return xmax_; }

  float get_ymax() 
  { return ymax_; }

  float get_zmax() 
  { return zmax_; }

  float get_lxmin() 
  { return lxmin_; }

  float get_lymin() 
  { return lymin_; }

  float get_lzmin() 
  { return lzmin_; }
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
  class_<Excitation>("Excitation", 
                     "Excitations applied to the FDTD grid", 
                     init<shared_ptr<Signal> >())
    //.def("excite", &ExcitationWrap::excite)
    //.def("excite", &ExcitationWrap::default_excite)
    .def("set_polarization", &Excitation::set_polarization)
    .def("set_type", &Excitation::set_type)
    .def("set_soft", &Excitation::set_soft)
    .def("get_soft", &Excitation::get_soft)
    .def("set_region", &Excitation::set_region)
    .def("set_time_param", &Excitation::set_time_param)
    ;

  class_<WindowedExcitation, WindowedExcitationWrap, bases<Excitation>,
    boost::noncopyable>("WindowedExcitation", 
                        "Excitations that apply a windowing function "
                        "to the excitation in the FDTD grid", 
                        init<shared_ptr<Signal> >())
    .def("init", &WindowedExcitation::init)
    .def("window", pure_virtual(&WindowedExcitation::window))

    .add_property("lxmin", &WindowedExcitationWrap::get_lxmin)
    .add_property("lymin", &WindowedExcitationWrap::get_lymin)
    .add_property("lzmin", &WindowedExcitationWrap::get_lzmin)

    .add_property("xmin", &WindowedExcitationWrap::get_xmin)
    .add_property("ymin", &WindowedExcitationWrap::get_ymin)
    .add_property("zmin", &WindowedExcitationWrap::get_zmin)

    .add_property("xmax", &WindowedExcitationWrap::get_xmax)
    .add_property("ymax", &WindowedExcitationWrap::get_ymax)
    .add_property("zmax", &WindowedExcitationWrap::get_zmax)
    //.def("excite", &WindowedExcitation::excite, 
    //     &WindowedExcitationWrap::default_excite)
    ;

  class_<Signal, SignalWrap, boost::noncopyable>
    ("Signal", "Make derived classes from this to create "
     "signal functions for excitations")
    .def("signal_function", &SignalWrap::signal_function)
    ;

  class_<Gaussm, bases<Signal> >("Gaussm", 
                                 "Gaussian modulated sine function")
    .def("signal_function", &Gaussm::signal_function)
    .def("set_parameters", &Gaussm::set_parameters)
    .def("get_alpha", &Gaussm::get_alpha)
    .def("get_deltaf", &Gaussm::get_deltaf)
    .def("get_f0", &Gaussm::get_f0)
    .def("length", &Gaussm::length)
    .def("signal_function", &Gaussm::signal_function) // in Signal
    ;

  class_<GaussPulse, bases<Signal> >("GaussPulse", 
                                     "Gaussian pulse function")
    .def("signal_function", &GaussPulse::signal_function)
    .def("set_parameters", &GaussPulse::set_parameters)
    .def("get_alpha", &GaussPulse::get_alpha)
    .def("get_taup", &GaussPulse::get_taup)
    .def("length", &GaussPulse::length)
    ;

  class_<DiffGaussPulse, bases<Signal> >("DiffGaussPulse", 
                                         "Differentiated Gaussian "
                                         "pulse function which peaks at 1/"
                                         "taup")
    .def("signal_function", &DiffGaussPulse::signal_function)
    .def("set_parameters", &DiffGaussPulse::set_parameters)
    .def("get_alpha", &DiffGaussPulse::get_alpha)
    .def("get_taup", &DiffGaussPulse::get_taup)
    .def("length", &DiffGaussPulse::length)
    ;

  class_<ExpSine, bases<Signal> >("ExpSine", "Ramping up sine function")
    .def(init<float>())
    .add_property("frequency", &ExpSine::get_frequency, 
                  &ExpSine::set_frequency)
    .add_property("amplitude", &ExpSine::get_amplitude, 
                  &ExpSine::set_amplitude)
    .def("signal_function", &ExpSine::signal_function)
    ;

  class_<BartlettExcitation, bases<WindowedExcitation> >
    ("BartlettExcitation", "Bartlett windowed excitation; an "
     "attempt at a plane wave.", 
     init<shared_ptr<Signal> >())
    .def("excite", &BartlettExcitation::excite)
    ;

  class_<WaveguideExcitation, bases<WindowedExcitation> >
    ("WaveguideExcitation", "Waveguide excitation; you know, "
     "for those pesky waveguides!", 
     init<shared_ptr<Signal> >())
    .def("excite", &WaveguideExcitation::excite)
    .def("set_mode", &WaveguideExcitation::set_mode)
    ;

  class_<GaussWindExcitation, bases<WindowedExcitation> >
    ("GaussWindow", "Gaussian windowed excitation; approximates "
     "a plane wave.", 
     init<shared_ptr<Signal> >())
    .def("excite", &GaussWindExcitation::excite)
    ;

  class_<PeriodicExcitation, bases<Excitation> >
    ("PeriodicExcitation", 
     "An excitation which can be used in conjunction with periodic boundary "
     "conditions.", init<shared_ptr<Signal> >())
    .def("excite", &PeriodicExcitation::excite)
    .def("set_region", (void(PeriodicExcitation::*)(shared_ptr<CSGBox>, Face))&PeriodicExcitation::set_region)
    .def("get_face", &PeriodicExcitation::get_face)
    ;
}
