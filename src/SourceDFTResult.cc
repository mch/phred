#include "SourceDFTResult.hh"
#include "Constants.hh"

SourceDFTResult::SourceDFTResult(TimeExcitation &te, 
                                  field_t freq_start,
                                  field_t freq_stop, 
                                  unsigned int num_freqs)
  : te_(te), freq_start_(freq_start), freq_stop_(freq_stop),
    num_freqs_(num_freqs)
{
  if (freq_stop_ < freq_start_)
  {
    field_t temp = freq_stop_;
    freq_stop_ = freq_start_;
    freq_start_ = temp;
  }

  time_dim_ = false; // We have only one output at the end. 

  dim_lens_.push_back(3); // Rows of (Freq, DFT)
  dim_lens_.push_back(num_freqs_); // num_freq rows
  var_name_ = "Source DFT";

  freq_space_ = (freq_stop_ - freq_start_) / num_freqs_;
  result_ = new field_t[(num_freqs_ + 1) * 3];

  for (unsigned int i = 0; i <= num_freqs_; i++)
  {
    result_[i * 3] = freq_start_ + freq_space_ * i;
    result_[i*3 + 1] = 0;
    result_[i*3 + 2] = 0;
  }
  
  MPI_Datatype temp;
  MPI_Type_contiguous(3, GRID_MPI_TYPE, &temp);
  MPI_Type_commit(&temp);

  data_.set_num(num_freqs_);
  data_.set_ptr(result_);
  data_.set_datatype(temp);
}

SourceDFTResult::~SourceDFTResult()
{ 
  delete[] result_;
}

void SourceDFTResult::set_excitation(const TimeExcitation &te)
{
  te_ = te;
}

Data &SourceDFTResult::get_result(Grid &grid, unsigned int time_step)
{
  field_t sf = te_.source_function(grid, time_step);
  for (unsigned int i = 0; i <= num_freqs_; i++)
  {
    result_[i*3 + 1] += sf * cos(2 * PI * result_[i*3] 
                                 * time_step * grid.get_deltat());

    result_[i*3 + 2] += (-1) * sf * sin(2 * PI * result_[i*3] 
                                        * time_step * grid.get_deltat());
  }

  return data_;
}
