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

#include "../DataWriter.hh"
#include "../AsciiDataWriter.hh"
#include "../Hdf5DataWriter.hh"
#include "../HdfDataWriter.hh"
#include "../NetCDFDataWriter.hh"
#include "../VtkDataWriter.hh"

using namespace boost::python;



BOOST_PYTHON_MODULE(datawriters)
{
  class_<NetCDFDataWriter, bases<DataWriter>, boost::noncopyable>("NetCDFDataWriter", "Writes result data to NetCDF files", init<int, int, const char *, bool>())
    .def(init<int, int, const char *>())
    .def(init<int, int>())
    .def("set_filename", &NetCDFDataWriter::set_filename)
    .def("add_variable", &NetCDFDataWriter::add_variable)
    ;

  class_<AsciiDataWriter, bases<DataWriter>, boost::noncopyable>("AsciiDataWriter", "Writes result data to ASCII files", init<int, int, const char *, Result &>())
    .def(init<int, int>())
    .def("set_filename", &AsciiDataWriter::set_filename)
    .def("add_variable", &AsciiDataWriter::add_variable)
    ;
}