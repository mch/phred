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

#include "../Results/Result.hh"
#include "../Results/PlaneResult.hh"
#include "../Results/PointResult.hh"
#include "../Results/PointDFTResult.hh"
#include "../Results/SignalDFTResult.hh"
#include "../Results/SignalTimeResult.hh"
#include "../Results/BlockResult.hh"
#include "../Results/FarfieldResult2.hh"
#include "../Results/FarfieldResult.hh"
#include "../Results/PowerResult.hh"
#include "../Results/AvgPowerResult.hh"
#include "../Results/FakeResult.hh"
#include "../Results/GridResult.hh"
#include "../Results/SurfaceCurrentResult.hh"
#include "../Results/WindowedExResult.hh"
#include "../Results/PropertiesResult.hh"

using namespace boost::python;

void export_results()
{
  class_<Result, boost::noncopyable>("Result", no_init)
    .def("set_time_param", &Result::set_time_param)
    .def("set_name", &Result::set_name)
    //.def("get_name", &Result::get_name)
    .def("set_dw_name", &Result::set_dw_name)
    ;

  class_<DFTResult, bases<Result>, boost::noncopyable>("DFTResult", no_init)
    .def("set_freq", &DFTResult::set_freq)
    .def("get_freq_start", &DFTResult::get_freq_start)
    .def("get_freq_stop", &DFTResult::get_freq_stop)
    .def("get_num_freq", &DFTResult::get_num_freq)
    ;

  class_<PlaneResult, bases<Result> >("PlaneResult")
    .def("set_plane", &PlaneResult::set_plane)
    .def("get_face", &PlaneResult::get_face)
    .def("set_field", &PlaneResult::set_field)
    ;

  class_<FakeResult, bases<Result> >("FakeResult")
    ;

  class_<GridResult, bases<Result> >("GridResult")
    ;

  class_<PointResult, bases<Result> >("PointResult")
    .def("set_point", &PointResult::set_point)
    .def("get_point", &PointResult::get_point)
    .add_property("point", &PointResult::get_point, &PointResult::set_point)
    ;

  class_<PointDFTResult, bases<DFTResult> >("PointDFTResult")
    .def("set_point", &PointDFTResult::set_point)
    .def("get_point", &PointDFTResult::get_point)
    .add_property("point", &PointDFTResult::get_point, 
                  &PointDFTResult::set_point)
    ;

  class_<PowerResult, bases<DFTResult> >("PowerResult")
    .def(init<field_t, field_t, unsigned int>())
    .def("set_region", &PowerResult::set_region)
    ;

  class_<AvgPowerResult, bases<DFTResult> >("AvgPowerResult")
    .def(init<field_t, field_t, unsigned int>())
    .def("set_region", &AvgPowerResult::set_region)
    ;

  class_<SignalDFTResult, bases<DFTResult> >
    ("SignalDFTResult", init<Signal &>())
    .def(init<Signal &, float, float, unsigned int>())
    ;

  class_<SignalTimeResult, bases<Result> >
    ("SignalTimeResult", init<Signal &>())
    ;

  class_<BlockResult, bases<Result> >
    ("BlockResult", "Outputs a single field component in part of the grid.")
    .def("set_region", &BlockResult::set_region)
    .def("get_region", &BlockResult::get_region)
    .def("set_field", &BlockResult::set_field)
    .def("get_field", &BlockResult::get_field)
    ;

//   class_<SurfaceCurrentResult, bases<Result> >
//     ("SurfaceCurrentResult", 
//      "Generates data about the currents excited on the surface "
//      "of an equivalent surface.")
//     .def("set_region", &SurfaceCurrentResult::set_region, 
//          "Set the CSGBox which defines the surface over which "
//          "currents should be calculated.")
//     ;

  class_<FarfieldResult2, bases<DFTResult> >("FarfieldResult2", 
                                             "Calculates the farfield "
                                             "resulting from a scattering "
                                             "object in the near field.")
    .def("set_theta", &FarfieldResult2::set_theta, 
         "Set the angle range from the X axis in radians")
    .def("set_phi", &FarfieldResult2::set_phi,
         "Set the angle range from the Z axis in radians")
    .def("set_theta_degrees", &FarfieldResult2::set_theta_degrees,
         "Set the angle range from the X axis in degrees")
    .def("set_phi_degrees", &FarfieldResult2::set_phi_degrees,
         "Set the angle range from the Z axis in degrees")
    .def("set_radius", &FarfieldResult2::set_radius, 
         "Set the radius of the sphere where observations are being done.")
    .def("get_radius", &FarfieldResult2::get_radius,
         "Get the radius of the sphere where observations are being done.")
    .def("set_region", &FarfieldResult2::set_region, 
         "Set the CSGBox which defines the surface over which "
         "currents should be calculated.")
    .def("use_face", &FarfieldResult2::use_face,
         "Exclude or include faces from the Huygen's surface. This breaks "
         "the surface equivalence theorm slightly, use with caution. ")
    ;

  class_<FarfieldResult, bases<DFTResult> >("FarfieldResult", 
                                            "Calculates the farfield "
                                            "resulting from a scattering "
                                            "object in the near field using "
                                            "Luebbers' time domain method.")
    .def("set_theta", &FarfieldResult::set_theta, 
         "Set the angle range from the X axis in radians")
    .def("set_phi", &FarfieldResult::set_phi,
         "Set the angle range from the Z axis in radians")
    .def("set_theta_degrees", &FarfieldResult::set_theta_degrees,
         "Set the angle range from the X axis in degrees")
    .def("set_phi_degrees", &FarfieldResult::set_phi_degrees,
         "Set the angle range from the Z axis in degrees")
    .def("set_radius", &FarfieldResult::set_radius, 
         "Set the radius of the sphere where observations are being done.")
    .def("get_radius", &FarfieldResult::get_radius,
         "Get the radius of the sphere where observations are being done.")
    .def("set_region", &FarfieldResult::set_region, 
         "Set the CSGBox which defines the surface over which "
         "currents should be calculated.")
    .def("use_face", &FarfieldResult::use_face,
         "Exclude or include faces from the Huygen's surface. This breaks "
         "the surface equivalence theorem slightly, use with caution. ")
    ;

//   class_<FarfieldResult, bases<Result> >("FarfieldResult")
//     .def("set_huygens", &FarfieldResult::set_huygens)
//     .def("set_output_type", &FarfieldResult::set_output_type)
//     .def("set_frequencies", &FarfieldResult::set_frequencies)
//     .def("set_angles", &FarfieldResult::set_angles)
//     .def("set_axis", &FarfieldResult::set_axis)
//     ;

  class_<WindowedExResult, bases<Result> >
    ("WindowedExResult", init<shared_ptr<WindowedExcitation> >())
    ;

  class_<PropertiesResult, bases<Result> >
    ("PropertiesResult", "Returns properties of the simulation such as "
     "the grid and time deltas.")
    ;
}
