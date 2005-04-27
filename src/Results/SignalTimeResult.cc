/* 
   phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004-2005 Matt Hughes <mhughe@uvic.ca>

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

#include "SignalTimeResult.hh"
#include "../Globals.hh"

SignalTimeResult::SignalTimeResult(Signal &te)
  : te_(te)
{
  variables_["Signal time excitation"] = &var_;
  var_.add_dimension("Signal Time Excitation", 2, 2, 0);

  MPI_Datatype temp;
  MPI_Type_contiguous(2, GRID_MPI_TYPE, &temp);
  MPI_Type_commit(&temp);

  var_.set_ptr(result_);
  var_.set_datatype(temp);
}

SignalTimeResult::~SignalTimeResult()
{ }

void SignalTimeResult::init(const Grid &grid)
{
  var_.set_name(base_name_);
}

void SignalTimeResult::set_excitation(const Signal &te)
{
  te_ = te;
}

void SignalTimeResult::calculate_result(const Grid &grid, 
                                        unsigned int time_step)
{
  if (result_time(time_step) && MPI_RANK == 0) 
  {
    result_[0] = grid.get_deltat() * time_step; 
    result_[1] = te_.signal_function(result_[0]);
    var_.set_num(1); 

//     cerr << "SignalTimeResult, data is " << result_[0] 
//          << " and " << result_[1] << endl;
//     cerr << "Pointer is " << reinterpret_cast<void *>(result_) << endl;

  } 
  else {
    var_.set_num(0); 
  }

}

ostream& SignalTimeResult::to_string(ostream &os) const
{
  return os << "SignalTimeResult";
}
