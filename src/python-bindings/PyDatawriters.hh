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
#include "../MatlabDataWriter.hh"
//#include "../Hdf5DataWriter.hh"
//#include "../HdfDataWriter.hh"
#include "../NetCDFDataWriter.hh"
//#include "../VtkDataWriter.hh"
//#include "../PyDataWriter.hh"

using namespace boost::python;

/**
 * Wrapper class for data writers
 */
class DataWriterWrap : public DataWriter
{
private:
  PyObject *self_;

public:
  DataWriterWrap(PyObject *self, int rank, int size)
    : self_(self), DataWriter(rank, size)
  {}

  unsigned int write_data(unsigned int time_step, 
                          Variable &var, MPI_Datatype t, 
                          void *ptr, unsigned int len)
  {}

  void add_variable(Result &result)
  {}

};

BOOST_PYTHON_MODULE(DataWriters)
{
  class_<DataWriter, DataWriterWrap, boost::noncopyable>("DataWriter", 
                                                         "Data writer base class",
                                                         init<int, int>())
    .def("set_filename", &DataWriter::set_filename)
    ;

  class_<NetCDFDataWriter, bases<DataWriter>, boost::noncopyable>("NetCDFDataWriter", 
                           "Writes result data to NetCDF files", 
                           init<int, int, const char *, bool>())
    .def(init<int, int, const char *>())
    .def(init<int, int>())
    .def("set_filename", &NetCDFDataWriter::set_filename)
    .def("add_variable", &NetCDFDataWriter::add_variable)
    ;

  class_<AsciiDataWriter, bases<DataWriter>, boost::noncopyable>("AsciiDataWriter", 
                          "Writes result data to ASCII files", 
                          init<int, int, const char *, Result &>())
    .def(init<int, int>())
    .def("set_filename", &AsciiDataWriter::set_filename)
    .def("add_variable", &AsciiDataWriter::add_variable)
    ;

  class_<MatlabDataWriter, bases<DataWriter>, boost::noncopyable>("MatlabDataWriter", 
                           "Writes data to Matlab 5 format files.",
                           init<int, int, const char *, Result &>())
    .def(init<int, int>())
    .def("add_variable", &MatlabDataWriter::add_variable)
    .def("test", &MatlabDataWriter::test)
    .def("set_filename", &MatlabDataWriter::set_filename)
    ;

//   class_<Hdf5DataWriter>("Hdf5DataWriter", 
//                          "Writes result data to HDF 5 files", 
//                           init<int, int, const char *, Result &>())
//     .def(init<int, int>())
//     .def("set_filename", &Hdf5DataWriter::set_filename)
//     .def("add_variable", &Hdf5DataWriter::add_variable)
//     ;

//   class_<HdfDataWriter>("HdfDataWriter", 
//                          "Writes result data to HDF files", 
//                           init<int, int, const char *, Result &>())
//     .def(init<int, int>())
//     .def("set_filename", &HdfDataWriter::set_filename)
//     .def("add_variable", &HdfDataWriter::add_variable)
//     ;

//   class_<VtkDataWriter>("VtkDataWriter", 
//                          "Writes result data to VTK files", 
//                           init<int, int, const char *, Result &>())
//     .def(init<int, int>())
//     .def("set_filename", &VtkDataWriter::set_filename)
//     .def("add_variable", &VtkDataWriter::add_variable)
//     ;

//   class_<PyDataWriter>("PyDataWriter", 
//                        "Exports data into Python so user scripts can work with it. "
//                        init<int, int>())
//     .def("add_variable", &PyDataWriter::add_variable)
//     ;
}
