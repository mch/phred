#include <boost/python.hpp>

#include "../FDTD.hh"

using namespace boost::python;


BOOST_PYTHON_MODULE(FDTD)
{
  class_<FDTD>("FDTD")
    .def("set_grid_size", &FDTD::set_grid_size)
    .def("set_grid_deltas", &FDTD::set_grid_deltas)
    .def("set_boundary", &FDTD::set_boundary)
    .def("load_materials", &FDTD::load_materials)
    .def("add_excitation", &FDTD::add_excitation)
    .def("add_result", &FDTD::add_result)
    .def("add_datawriter", &FDTD::add_datawriter)
    .def("map_result_to_dw", &FDTD::map_result_to_datawriter)
    .def("run", &FDTD::run)
    ;

  
}
