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
field_t call_sf(SourceFunction& sf, Grid &grid, 
                    unsigned int time_step) 
{ return sf.source_function(grid, time_step); }

/**
 * Since the normal C++ version of Excitation is a template class,
 * this is instance of that template which stores a reference to
 * SourceFunction, the base class for these things. A proxy class is
 * exposed which takes a pointer to a subclass of SourceFunction so
 * that it can be used from Python. This wrapper also allows for
 * derived classes built in Python. 
 */
class ExcitationWrap : public Excitation<SourceFunction>
{
  PyObject* self_;

public:
  ExcitationWrap(PyObject* self, SourceFunction &sf)
    :  Excitation<SourceFunction>(sf), self_(self) {}

  void excite(Grid &grid, unsigned int time_step,
              FieldType type) 
  { return call_method<void>(self, "excite"); }

  void default_excite(Grid &grid, unsigned int time_step,
                      FieldType type) 
  { Excitation::excite(grid, time_step,
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

BOOST_PYTHON_MODULE(excitation)
{
  class_<Excitation, ExcitationWrap, boost::noncopyable>("Excitation")
    .def("call_excite", &Excitation::excite, 
         &ExcitationWrap::default_excite)
    ;

  class_<SourceFunction, SourceFunctionWrap, boost::noncopyable>("SourceFunction", no_init)
    .def("call_sf", call_sf);

  class_<Gaussm, bases<SourceFunction> >("Gaussm")
    .def("set_parameters", &Gaussm::set_parameters)
    .def("get_alpha", &Gaussm::get_alpha)
    .def("get_deltaf", &Gaussm::get_deltaf)
    .def("get_f0", &Gaussm::get_f0)
    .def("source_function", call_sf)
    .def("excite", call_excite)
    ;
}
