#include <boost/python.hpp>

#include "../Grid.hh"

using namespace boost::python;

class GridWrap : public Grid
{
  PyObject* self;

  GridWrap(PyObject* self_)
    : self(self_) {}

  void apply_boundaries() 
  { 
    return call_apply_boundaries<void>(self, "apply_boundaries"); 
  }

  void load_materials() 
  { 
    return call_load_materials<void>(self, "load_materials"); 
  }

  void free_material() 
  { 
    return call_free_material<void>(self, "free_material"); 
  }

  void setup_grid() 
  { 
    return call_setup_grid<void>(self, "setup_grid"); 
  }

  void alloc_grid() 
  { 
    return call_alloc_grid<void>(self, "alloc_grid"); 
  }

  void free_grid() 
  { 
    return call_free_grid<void>(self, "free_grid"); 
  }

  void update_ex() 
  { return call_update_ex<void>(self, "update_ex"); }

  void update_ey() 
  { return call_update_ey<void>(self, "update_ey"); }

  void update_ez() 
  { return call_update_ez<void>(self, "update_ez"); }

  void update_hx() 
  { return call_update_hx<void>(self, "update_hx"); }

  void update_hy() 
  { return call_update_hy<void>(self, "update_hy"); }

  void update_hz() 
  { return call_update_hz<void>(self, "update_hz"); }
};

BOOST_PYTHON_MODULE(Grid)
{
  class_<Grid>("Grid")
    .def(init<const Grid &>())
    .def("greet", &World::greet)
    .def("set", &World::set)
    ;

  class_<Grid, GridWrap, boost::noncopyable>("Grid", no_init)
    .def("call_f", call_f);
}
