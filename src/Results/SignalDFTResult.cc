/* 
   phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004 Matt Hughes <mhughe@uvic.ca>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#include "SignalDFTResult.hh"
#include "../Constants.hh"
#include "../Globals.hh"

#include <cmath>

SignalDFTResult::SignalDFTResult(Signal &te)
  : te_(te), result_(0)
{
  variables_["Signal_DFT"] = &var_;
}

SignalDFTResult::SignalDFTResult(Signal &te, 
                                  field_t freq_start,
                                  field_t freq_stop, 
                                  unsigned int num_freqs)
  : DFTResult(freq_start, freq_stop, num_freqs), te_(te), 
    result_(0)
{
  variables_["Signal_DFT"] = &var_;
}

SignalDFTResult::~SignalDFTResult()
{ 
  if (result_)
    delete[] result_;
}

void SignalDFTResult::init(const Grid &grid)
{
  if (frequencies_.length() == 0)
  {
    throw ResultException("SignalDFTResult: Frequency interval is not"
                          " correctly set up. ");
  }

  var_.has_time_dimension(false); // We have only one output at the end. 
  var_.add_dimension("freqs", frequencies_.length(), frequencies_.length(), 0);
  var_.add_dimension("results", 3, 3, 0);
  var_.set_name(base_name_);

  result_ = new field_t[(frequencies_.length() + 1) * 3];

  if (!result_)
    throw MemoryException();
  
  for (unsigned int i = 0; i <= frequencies_.length(); i++)
  {
    result_[i * 3] = frequencies_.get(i); // _start() + freq_space_ * i;
    result_[i*3 + 1] = 0;
    result_[i*3 + 2] = 0;
  }
  
  MPI_Datatype temp;
  MPI_Type_contiguous(3, GRID_MPI_TYPE, &temp);
  MPI_Type_commit(&temp);

  var_.set_num(frequencies_.length());
  var_.set_ptr(result_);
  var_.set_datatype(temp);
}

void SignalDFTResult::deinit()
{
  if (result_) {
    delete[] result_;
    result_ = 0;
  }
}

void SignalDFTResult::set_excitation(const Signal &te)
{
  te_ = te;
}

void SignalDFTResult::calculate_result(const Grid &grid, 
                                       unsigned int time_step)
{
  if (result_time(time_step) && MPI_RANK == 0) 
  {
    float time = time_step * grid.get_deltat();
    field_t sf = te_.signal_function(time);
    for (unsigned int i = 0; i <= frequencies_.length(); i++)
    {
      result_[i*3 + 1] += sf * cos(-2 * PI * result_[i*3] 
                                   * time);
    
      result_[i*3 + 2] += (-1) * sf * sin(-2 * PI * result_[i*3] 
                                          * time);
    }
    var_.set_num(frequencies_.length()); 
  } else {
    var_.set_num(0); 
  }

}

ostream& SignalDFTResult::to_string(ostream &os) const
{
  os << "SignalDFTResult outputting " << frequencies_.length() 
     << ", starting at " << frequencies_.get_start() << " and ending at "
     << frequencies_.get_end();

  return os;
}
