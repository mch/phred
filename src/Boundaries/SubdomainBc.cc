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


#include "SubdomainBc.hh"
#include "Grid.hh"

void SubdomainBc::apply(Face face, Grid &grid, FieldType type)
{

  // This has been changed to use the RxTx interface, set up in Grid.
//   region_t r = find_face(face, grid);
//   MPI_Datatype t;

//   t = grid.get_plane_dt(face);

//   // Send away!
//   if (type == E) {
//     send_recv(grid.get_face_start(face, FC_EX, 1),
//               grid.get_face_start(face, FC_EX, 0), t);
//     send_recv(grid.get_face_start(face, FC_EY, 1),
//               grid.get_face_start(face, FC_EY, 0), t);
//     send_recv(grid.get_face_start(face, FC_EZ, 1),
//               grid.get_face_start(face, FC_EZ, 0), t);
//   } else if (type == H) {
//     send_recv(grid.get_face_start(face, FC_HX, 1),
//               grid.get_face_start(face, FC_HX, 0), t);
//     send_recv(grid.get_face_start(face, FC_HY, 1),
//               grid.get_face_start(face, FC_HY, 0), t);
//     send_recv(grid.get_face_start(face, FC_HZ, 1),
//               grid.get_face_start(face, FC_HZ, 0), t);
//   }

  // Send and recieve and data we've be contracted to do. 
  vector<RxTxData>::iterator iter = rx_tx_data_.begin();
  vector<RxTxData>::iterator iter_end = rx_tx_data_.end();

  while (iter != iter_end) 
  {
    if (type == (*iter).get_field_type()) 
    {
      MPI_Datatype dtype = (*iter).get_datatype();
      send_recv((*iter).get_tx_ptr(), (*iter).get_rx_ptr(), 
                dtype);
    }

    ++iter;
  }
}

void SubdomainBc::send_recv(void *tx_ptr, void *rx_ptr, MPI_Datatype &t)
{
  MPI_Status status;

  if (rank_ > neighbour_) {
    MPI_Send(tx_ptr, 1, t, neighbour_, 1, MPI_COMM_WORLD);
    MPI_Recv(rx_ptr, 1, t, neighbour_, 1, MPI_COMM_WORLD, &status);
  } else {
    MPI_Recv(rx_ptr, 1, t, neighbour_, 1, MPI_COMM_WORLD, &status);
    MPI_Send(tx_ptr, 1, t, neighbour_, 1, MPI_COMM_WORLD);
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
