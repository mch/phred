#include <boost/python.hpp>

#include "../Result.hh"
#include "../PlaneResult.hh"
#include "../PointResult.hh"
#include "../PointDFTResult.hh"
#include "../SourceDFTResult.hh"
#include "../SourceTimeResult.hh"
#include "../BlockResult.hh"

using namespace boost::python;

/**
 * Free function to facilitate ResultWrap; calling get_result on a
 * Python object derived from Result.
 */
Data &call_get_result(Result &r, const Grid &grid, unsigned int time_step)
{ return r.get_result(grid, time_step); }

/**
 * This allows for results written in Python
 */
class ResultWrap : public Result
{
private:
  PyObject *self_;

public:
  ResultWrap(PyObject *self)
    : self_(self) {}

  Data &get_result(const Grid &grid, unsigned int time_step)
  { return call_method<Data &>(self, "get_result"); }

};

BOOST_PYTHON_MODULE(results)
{
  class_<Result, ResultWrap, boost::noncopyable>("Result")
    .def("has_time_dimension", &Result::has_time_dimension)
    .def("set_time_param", &Result::set_time_param)
    .def("set_dw_name", &Result::set_dw_name)
    .def("get_dw_name", &Result::get_dw_name)
    .def("set_name", &Result::set_name)
    .def("get_name", &Result::get_name)
    .def("get_dim_lengths", &Result::get_dim_lengths)
    .def("get_dim_names", &Result::get_dim_names)
    ;

  def("call_get_result", call_get_result);

  class_<PlaneResult, bases<Result> >("PlaneResult")
    .def("set_plane", &PlaneResult::set_plane)
    .def("get_plane", &PlaneResult::get_plane)
    .def("get_face", &PlaneResult::get_face)
    .def("set_face", &PlaneResult::set_face)
    .def("set_field", &PlaneResult::set_field)
    ;

  class_<PointResult, bases<Result> >("PointResult")
    .def("set_point", &PointResult::set_point)
    .def("get_point", &PointResult::get_point)
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
    ;

  class_<SourceDFTResult, bases<Result> >("SourceDFTResult")
    .def("set_freq_start", &SourceDFTResult::set_freq_start)
    .def("get_freq_start", &SourceDFTResult::get_freq_start)
    .def("set_freq_stop", &SourceDFTResult::set_freq_stop)
    .def("get_freq_stop", &SourceDFTResult::get_freq_stop)
    .def("set_num_freq", &SourceDFTResult::set_num_freq)
    .def("get_num_freq", &SourceDFTResult::get_num_freq)
    ;

  class_<SourceTimeResult, bases<Result> >("SourceTimeResult")
    ;

  class_<BlockResult, bases<Result> >("BlockResult")
    .def("set_region", &BlockResult::set_region)
    .def("get_region", &BlockResult::get_region)
    .def("set_field_component", &BlockResult::set_field_component)
    .def("get_field_component", &BlockResult::get_field_component)
    ;
}
