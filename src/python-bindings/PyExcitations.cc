#include <boost/python.hpp>

#include "../Gaussm.hh"

using namespace boost::python;

class TimeExcitationWrap : public TimeExcitation
{
  PyObject* self;

  TimeExcitationWrap(PyObject* self_)
    : self(self_) {}

  void excite(Grid &grid, unsigned int time_step,
              FieldType type) 
  { return call_excite<void>(self, "excite"); }

  field_t source_function(Grid &grid, unsigned int time_step) 
  { return call_source_function<field_t>(self, "source_function"); }

  int default_excite(Grid &grid, unsigned int time_step,
                     FieldType type) 
  { return TimeExcitation::excite(grid, time_step,
                                  type); }
};

BOOST_PYTHON_MODULE(Excitation)
{
  class_<TimeExcitation, TimeExcitationWrap, boost::noncopyable>("TimeExcitation", no_init)
    .def("call_excite", &TimeExcitation::excite, 
         &TimeExcitationWrap::default_excite)
    .def("call_source_function", call_source_function)
    ;
}
