#include <boost/python.hpp>

#include "../Gaussm.hh"

using namespace boost::python;

/**
 * Since the normal C++ version of Excitation is a template class,
 * this is instance of that template which stores a reference to
 * SourceFunction, the base class for these things. A proxy class is
 * exposed which takes a pointer to a subclass of SourceFunction so
 * that it can be used from Python. 
 */
class PyExcitation : public Excitation<SourceFunction> 
{
protected:
public:
  PyExcitation(SourceFunction &sf)
    : Excitation<SourceFunction>(sf)
  {}
  ~PyExcitation();
};

void call_excite(TimeExcitation& te, Grid &grid, 
                        unsigned int time_step, FieldType type) 
{ return te.excite(grid, time_step, type); }

field_t call_sf(TimeExcitation& te, Grid &grid, 
                    unsigned int time_step) 
{ return te.source_function(grid, time_step); }


class TimeExcitationWrap : public TimeExcitation
{
  PyObject* self;

public:
  TimeExcitationWrap(PyObject* self_)
    : self(self_) {}

  void excite(Grid &grid, unsigned int time_step,
              FieldType type) 
  { return call_method<void>(self, "excite"); }

  field_t source_function(Grid &grid, unsigned int time_step) 
  { return call_method<field_t>(self, "source_function"); }

  void default_excite(Grid &grid, unsigned int time_step,
                     FieldType type) 
  { TimeExcitation::excite(grid, time_step,
                           type); }
};

BOOST_PYTHON_MODULE(excitation)
{
  class_<TimeExcitation, TimeExcitationWrap, boost::noncopyable>("TimeExcitation", no_init)
    .def("call_excite", &TimeExcitation::excite, 
         &TimeExcitationWrap::default_excite)
    .def("call_sf", call_sf)
    ;

  class_<Gaussm, bases<TimeExcitation> >("Gaussm")
    .def("set_parameters", &Gaussm::set_parameters)
    .def("get_alpha", &Gaussm::get_alpha)
    .def("get_deltaf", &Gaussm::get_deltaf)
    .def("get_f0", &Gaussm::get_f0)
    .def("source_function", call_sf)
    .def("excite", call_excite)
    ;
}
