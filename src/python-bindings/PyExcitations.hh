#include <boost/python.hpp>

#include "../Excitation.hh"
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
field_t call_source_functionf(SourceFunction& sf, Grid &grid, 
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

  void excite(Grid &grid, unsigned int time_step,
              FieldType type) 
  { return call_method<void>(self, "excite"); }

  void default_excite(Grid &grid, unsigned int time_step,
                      FieldType type) 
  { Excitation::excite(grid, time_step,
                       type); }
};

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
  { return call_method<void>(self, "excite"); }

  void default_excite(Grid &grid, unsigned int time_step,
                      FieldType type) 
  { WindowedExcitation::excite(grid, time_step,
                               type); }  
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

  field_t source_function(Grid &grid, unsigned int time_step)
  { return call_method<field_t>(self, "source_function"); }
};

BOOST_PYTHON_MODULE(excitations)
{
  class_<Excitation, ExcitationWrap, boost::noncopyable>("Excitation")
    .def("excite", &Excitation::excite, 
         &ExcitationWrap::default_excite)
    .def("set_polarization", &Excitation::set_polatization)
    .def("set_type", &Excitation::set_type)
    .def("set_soft", &Excitation::set_soft)
    .def("get_soft", &Excitation::get_soft)
    .def("set_region", &Excitation::set_region)
    ;

  class_<WindowedExcitation, WindowedExcitationWrap, bases<Excitation>, boost:noncopyable>("WindowedExcitation")
    .def("excite", &WindowedExcitation::excite, 
         &WindowedExcitationWrap::default_excite)
    ;

  class_<SourceFunction, SourceFunctionWrap, boost::noncopyable>("SourceFunction")
    .def("source_function", &SourceFunction::source_function)
    ;
  //    .def("call_sf", call_sf)

  class_<Gaussm, bases<SourceFunction> >("Gaussm")
    .def("set_parameters", &Gaussm::set_parameters)
    .def("get_alpha", &Gaussm::get_alpha)
    .def("get_deltaf", &Gaussm::get_deltaf)
    .def("get_f0", &Gaussm::get_f0)
    .def("source_function", &Gaussm::source_function)
    ;
}
