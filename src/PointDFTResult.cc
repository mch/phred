#include "PointDFTResult.hh"
#include "Constants.hh"

PointDFTResult::PointDFTResult(field_t freq_start,
                               field_t freq_stop, 
                               unsigned int num_freqs)
  : freq_start_(freq_start), freq_stop_(freq_stop),
    num_freqs_(num_freqs)
{
  if (freq_stop_ < freq_start_)
  {
    field_t temp = freq_stop_;
    freq_stop_ = freq_start_;
    freq_start_ = temp;
  }

  time_dim_ = false; // We have only one output at the end. 

  dim_lens_.push_back(13); // Rows of 1 col for freq, 2 cols each component
  dim_lens_.push_back(num_freqs_); // num_freq rows
  var_name_ = "Point DFT";

  freq_space_ = (freq_stop_ - freq_start_) / num_freqs_;
  result_ = new field_t[(num_freqs_ + 1) * 13];

  for (unsigned int i = 0; i <= num_freqs_; i++)
  {
    result_[i * 13] = freq_start_ + freq_space_ * i;
    for (unsigned int j = 1; j < 13; j++)
      result_[i*13 + j] = 0;
  }
  
  MPI_Datatype temp;
  MPI_Type_contiguous(13, GRID_MPI_TYPE, &temp);
  MPI_Type_commit(&temp);

  data_.set_num(num_freqs_);
  data_.set_ptr(result_);
  data_.set_datatype(temp);
}

PointDFTResult::~PointDFTResult()
{ 
  delete[] result_;
}

Data &PointDFTResult::get_result(Grid &grid, unsigned int time_step)
{
  delta_t dt = grid.get_deltat();
  delta_t time = dt * time_step;

  for (unsigned int i = 0; i <= num_freqs_; i++)
  {
    result_[i*13 + 1] += grid.get_ex(point_.x, point_.y, point_.z)
      * cos(2 * PI * result_[i*13] * time);
    
    result_[i*13 + 2] += (-1) * grid.get_ex(point_.x, point_.y, point_.z) 
      * sin(2 * PI * result_[i*13] * time);

    result_[i*13 + 3] += grid.get_ey(point_.x, point_.y, point_.z)
      * cos(2 * PI * result_[i*13] * time);
    
    result_[i*13 + 4] += (-1) * grid.get_ey(point_.x, point_.y, point_.z) 
      * sin(2 * PI * result_[i*13] * time);

    result_[i*13 + 5] += grid.get_ez(point_.x, point_.y, point_.z)
      * cos(2 * PI * result_[i*13] * time);
    
    result_[i*13 + 6] += (-1) * grid.get_ez(point_.x, point_.y, point_.z) 
      * sin(2 * PI * result_[i*13] * time);

    // H components
    result_[i*13 + 7] += grid.get_hx(point_.x, point_.y, point_.z)
      * cos(2 * PI * result_[i*13] * time);
    
    result_[i*13 + 8] += (-1) * grid.get_hx(point_.x, point_.y, point_.z) 
      * sin(2 * PI * result_[i*13] * time);

    result_[i*13 + 9] += grid.get_hy(point_.x, point_.y, point_.z)
      * cos(2 * PI * result_[i*13] * time);
    
    result_[i*13 + 10] += (-1) * grid.get_hy(point_.x, point_.y, point_.z) 
      * sin(2 * PI * result_[i*13] * time);

    result_[i*13 + 11] += grid.get_hz(point_.x, point_.y, point_.z)
      * cos(2 * PI * result_[i*13] * time);
    
    result_[i*13 + 12] += (-1) * grid.get_hz(point_.x, point_.y, point_.z) 
      * sin(2 * PI * result_[i*13] * time);
  }

  return data_;
}
