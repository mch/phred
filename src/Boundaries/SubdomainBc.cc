/* 
   Phred - Phred is a parallel finite difference time domain
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


#include "SubdomainBc.hh"
#include "../Grid.hh"
#include "../Globals.hh"

void SubdomainBc::apply(Face face, Grid &grid, FieldType type)
{
  // Send and recieve and data we've be contracted to do. 
  vector<RxTxData>::iterator iter = rx_tx_data_.begin();
  vector<RxTxData>::iterator iter_end = rx_tx_data_.end();

  while (iter != iter_end) 
  {
    if (type == (*iter).get_field_type()) 
    {
      MPI_Datatype dtype = (*iter).get_datatype();
      send_recv((*iter).get_tx_ptr(), 
                (*iter).get_rx_ptr(), 
                dtype);

// #ifdef DEBUG
//       int sz;
//       MPI_Aint ext;
//       MPI_Type_size((*iter).get_datatype(), &sz);
//       MPI_Type_extent((*iter).get_datatype(), &ext);
      
//       cout << "SubdomainBc::apply(), extent: " << ext << ", sent " 
//            << sz << " bytes from " 
//            << (*iter).get_tx_ptr() << "\nrecieved same number at " 
//            << (*iter).get_rx_ptr() << endl;
// #endif
    }

    ++iter;
  }
}

void SubdomainBc::send_recv(const void *tx_ptr, void *rx_ptr, MPI_Datatype &t)
{
  MPI_Status status;

  if (rank_ > neighbour_) {
    MPI_Send(const_cast<void *>(tx_ptr), 1, t, neighbour_, 1, MPI_COMM_PHRED);
    MPI_Recv(rx_ptr, 1, t, neighbour_, 1, MPI_COMM_PHRED, &status);
  } else {
    MPI_Recv(rx_ptr, 1, t, neighbour_, 1, MPI_COMM_PHRED, &status);
    MPI_Send(const_cast<void *>(tx_ptr), 1, t, neighbour_, 1, MPI_COMM_PHRED);
  }
}

void SubdomainBc::set_neighbour(int neighbour)
{
  neighbour_ = neighbour;
}

int SubdomainBc::get_neighbour()
{
  return neighbour_;
}

void SubdomainBc::set_rank(int rank)
{
  rank_ = rank;
}

int SubdomainBc::get_rank()
{
  return rank_;
}

void SubdomainBc::add_tx_rx_data(const RxTxData &x)
{
  rx_tx_data_.push_back(x);
}

BoundaryCondition SubdomainBc::get_type() const
{
  return SUBDOMAIN;
}
