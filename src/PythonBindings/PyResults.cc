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

#include "../Results/Result.hh"
#include "../Results/PlaneResult.hh"
#include "../Results/PointResult.hh"
#include "../Results/PointDFTResult.hh"
#include "../Results/SourceDFTResult.hh"
#include "../Results/SourceTimeResult.hh"
#include "../Results/BlockResult.hh"
#include "../Results/FarfieldResult.hh"
#include "../Results/PowerResult.hh"
#include "../Results/FakeResult.hh"

using namespace boost::python;


void export_results()
{
  class_<Result, boost::noncopyable>("Result", no_init)
    .def("set_time_param", &Result::set_time_param)
    .def("set_name", &Result::set_name)
    //.def("get_name", &Result::get_name)
    .def("set_dw_name", &Result::set_dw_name)
    ;

  class_<PlaneResult, bases<Result> >("PlaneResult")
    .def("set_plane", &PlaneResult::set_plane)
    .def("get_plane", &PlaneResult::get_plane)
    .def("get_face", &PlaneResult::get_face)
    .def("set_field", &PlaneResult::set_field)
    ;

  class_<FakeResult, bases<Result> >("FakeResult")
    ;

  class_<PointResult, bases<Result> >("PointResult")
    .def("set_point", &PointResult::set_point)
    .def("get_point", &PointResult::get_point)
    .add_property("point", &PointResult::get_point, &PointResult::set_point)
    ;

  class_<PointDFTResult, bases<Result> >("PointDFTResult")
    .def("set_point", &PointDFTResult::set_point)
    .def("get_point", &PointDFTResult::get_point)
    .def("set_freq_start", &PointDFTResult::set_freq_start)
    .def("get_freq_start", &PointDFTResult::get_freq_start)
    .def("set_freq_stop", &PointDFTResult::set_freq_stop)
    .def("get_freq_stop", &PointDFTResult::get_freq_stop)
    .def("set_num_freq", &PointDFTResult::set_num_freq)
    .def("get_num_freq", &PointDFTResult::get_num_freq)
    .add_property("freq_start", &PointDFTResult::get_freq_start,
                  &PointDFTResult::set_freq_start)
    .add_property("freq_stop", &PointDFTResult::get_freq_stop,
                  &PointDFTResult::set_freq_stop)
    .add_property("num_freqs", &PointDFTResult::get_num_freq,
                  &PointDFTResult::set_num_freq)
    .add_property("point", &PointDFTResult::get_point, 
                  &PointDFTResult::set_point)
    ;

  class_<PowerResult, bases<Result> >("PowerResult")
    .def(init<field_t, field_t, unsigned int>())
    .def("set_region", &PowerResult::set_region)
    .def("set_freq_start", &PowerResult::set_freq_start)
    .def("get_freq_start", &PowerResult::get_freq_start)
    .def("set_freq_stop", &PowerResult::set_freq_stop)
    .def("get_freq_stop", &PowerResult::get_freq_stop)
    .def("set_num_freq", &PowerResult::set_num_freq)
    .def("get_num_freq", &PowerResult::get_num_freq)
    .add_property("freq_start", &PowerResult::get_freq_start,
                  &PowerResult::set_freq_start)
    .add_property("freq_stop", &PowerResult::get_freq_stop,
                  &PowerResult::set_freq_stop)
    .add_property("num_freqs", &PowerResult::get_num_freq,
                  &PowerResult::set_num_freq)
    ;

  class_<SourceDFTResult, bases<Result> >("SourceDFTResult", init<SourceFunction &>())
    .def(init<SourceFunction &, float, float, unsigned int>())
    .def("set_freq_start", &SourceDFTResult::set_freq_start)
    .def("get_freq_start", &SourceDFTResult::get_freq_start)
    .def("set_freq_stop", &SourceDFTResult::set_freq_stop)
    .def("get_freq_stop", &SourceDFTResult::get_freq_stop)
    .def("set_num_freq", &SourceDFTResult::set_num_freq)
    .def("get_num_freq", &SourceDFTResult::get_num_freq)
    .add_property("freq_start", &SourceDFTResult::get_freq_start,
                  &SourceDFTResult::set_freq_start)
    .add_property("freq_stop", &SourceDFTResult::get_freq_stop,
                  &SourceDFTResult::set_freq_stop)
    .add_property("num_freqs", &SourceDFTResult::get_num_freq,
                  &SourceDFTResult::set_num_freq)
    ;

  class_<SourceTimeResult, bases<Result> >("SourceTimeResult", init<SourceFunction &>())
    ;

  class_<BlockResult, bases<Result> >("BlockResult")
    .def("set_region", &BlockResult::set_region)
    .def("get_region", &BlockResult::get_region)
    .def("set_field", &BlockResult::set_field)
    .def("get_field", &BlockResult::get_field)
    ;

  class_<FarfieldResult, bases<Result> >("FarfieldResult")
    .def("set_huygens", &FarfieldResult::set_huygens)
    .def("set_output_type", &FarfieldResult::set_output_type)
    .def("set_frequencies", &FarfieldResult::set_frequencies)
    .def("set_angles", &FarfieldResult::set_angles)
    .def("set_axis", &FarfieldResult::set_axis)
    ;
}
