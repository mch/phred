#include <boost/python.hpp>

#include "../Gaussm.hh"

using namespace boost::python;

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

BOOST_PYTHON_MODULE(Excitation)
{
  class_<TimeExcitation, TimeExcitationWrap, boost::noncopyable>("TimeExcitation", no_init)
    .def("call_excite", &TimeExcitation::excite, 
         &TimeExcitationWrap::default_excite)
    .def("call_sf", call_sf)
    ;
}
