#include "PlaneResult.hh"

PlaneResult::PlaneResult()
  : face_(FRONT), field_(EY)
{
  plane_.x = 0;
  plane_.y = 0;
  plane_.z = 0;

  num_dims_ = 2;
  dim_lens_ = new unsigned int[2];
  dim_lens_[0] = 0;
  dim_lens_[1] = 0;
}

PlaneResult::~PlaneResult()
{
  delete[] dim_lens_;
}

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

void PlaneResult::set_size(unsigned int x, unsigned int y)
{
  dim_lens_[0] = x;
  dim_lens_[1] = y;
}
