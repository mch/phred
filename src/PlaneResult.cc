#include "PlaneResult.hh"

PlaneResult::PlaneResult()
  : face_(FRONT), field_(EY)
{
  plane_.x = 0;
  plane_.y = 0;
  plane_.z = 0;
}

PlaneResult::~PlaneResult()
{}

Data &PlaneResult::get_result(Grid &grid, unsigned int time_step)
{
  if (result_time(time_step))
  {
    data_.set_num(1);
    data_.set_ptr(grid.get_face_start(face_, field_, plane_));

    MPI_Datatype t = grid.get_plane_dt(face_);
    data_.set_datatype(t);
  } else {
    data_.set_num(0);
  }

  return data_;
}
