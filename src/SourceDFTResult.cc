#include "SourceDFTResult.hh"

SourceDFTResult::SourceDFTResult(TimeExcitation &te, 
                                  field_t freq_start,
                                  field_t freq_stop, 
                                  unsigned int num_freqs)
  : te_(te), freq_start_(freq_start), freq_stop_(freq_stop),
    num_freqs_(num_freqs)
{
  dim_lens_.push_back(num_freqs);
  var_name_ = "Source DFT";

  MPI_Datatype temp;
  MPI_Type_contiguous(2, GRID_MPI_TYPE, &temp);
  MPI_Type_commit(&temp);

  data_.set_num(1);
  data_.set_ptr(result_);
  data_.set_datatype(temp);
}

SourceDFTResult::~SourceDFTResult()
{ }

void SourceDFTResult::set_excitation(const TimeExcitation &te)
{
  te_ = te;
}

Data &SourceDFTResult::get_result(Grid &grid, unsigned int time_step)
{
  result_[0] = grid.get_deltat() * time_step; 
  result_[1] = te_.source_function(grid, time_step);

  return data_;
}
