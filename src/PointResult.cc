#include "PointResult.hh"

PointResult::PointResult()
{
  MPI_Datatype temp;
  MPI_Type_contiguous(7, GRID_MPI_TYPE, &temp);
  MPI_Type_commit(&temp);
  data_.set_datatype(temp);

  num_dims_ = 1;
  dim_lens_ = new unsigned int[1];
  dim_lens_[0] = 7;
}

PointResult::~PointResult()
{
  delete[] dim_lens_;
}

Data &PointResult::get_result(Grid &grid, unsigned int time_step)
{
  point_t l;
  bool ours = true;

  if (point_.x < grid.get_lsx() 
      || point_.x >= grid.get_lsx() + grid.get_ldx())
    ours = false;
  else 
    l.x = point_.x - grid.get_lsx();

  if (point_.y < grid.get_lsy() 
      || point_.y >= grid.get_lsy() + grid.get_ldy())
    ours = false;
  else 
    l.y = point_.y - grid.get_lsy();

  if (point_.z < grid.get_lsz() 
      || point_.z >= grid.get_lsz() + grid.get_ldz())
    ours = false;
  else 
    l.z = point_.z - grid.get_lsz();

  if (ours) 
  {
    field_data_[0] = grid.get_deltat() * time_step;
    field_data_[1] = grid.get_ex(l.x, l.y, l.z);
    field_data_[2] = grid.get_ey(l.x, l.y, l.z);
    field_data_[3] = grid.get_ez(l.x, l.y, l.z);
    
    field_data_[4] = grid.get_hx(l.x, l.y, l.z);
    field_data_[5] = grid.get_hy(l.x, l.y, l.z);
    field_data_[6] = grid.get_hz(l.x, l.y, l.z);
    
    data_.set_ptr(field_data_);
    data_.set_num(1);
  } 
  else
  {
    data_.set_ptr(0);
    data_.set_num(0);
  }
  
  return data_;
}
