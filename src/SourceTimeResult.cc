#include "SourceTimeResult.hh"

SourceTimeResult::SourceTimeResult(SourceFunction &te)
  : te_(te)
{
  dim_lens_.push_back(2);
  var_name_ = "Source Time Excitation";

  MPI_Datatype temp;
  MPI_Type_contiguous(2, GRID_MPI_TYPE, &temp);
  MPI_Type_commit(&temp);

  data_.set_ptr(result_);
  data_.set_datatype(temp);
}

SourceTimeResult::~SourceTimeResult()
{ }

void SourceTimeResult::set_excitation(const SourceFunction &te)
{
  te_ = te;
}

Data &SourceTimeResult::get_result(const Grid &grid, unsigned int time_step)
{
  if (result_time(time_step)) 
  {
    result_[0] = grid.get_deltat() * time_step; 
    result_[1] = te_.source_function(grid, time_step);
    data_.set_num(1); 
  } 
  else {
    data_.set_num(0); 
  }

  return data_;
}
