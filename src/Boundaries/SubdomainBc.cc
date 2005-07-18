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

SubdomainBc::SubdomainBc()
  : locked_(false), num_(0), num_used_(0), req_(0)
{}

SubdomainBc::~SubdomainBc()
{
  if (req_)
  {
    delete[] req_;
    req_ = 0;
  }
}

void SubdomainBc::sd_deinit()
{
  if (req_)
  {
    delete[] req_;
    req_ = 0;
  }
}

void SubdomainBc::sd_init(const Grid &grid, Face face)
{}

void SubdomainBc::apply(Face face, Grid &grid, FieldType type)
{
  // Send and recieve and data we've be contracted to do. 
  vector<RxTxData>::iterator iter = rx_tx_data_.begin();
  vector<RxTxData>::iterator iter_end = rx_tx_data_.end();
  int idx = 0;

  while (iter != iter_end) 
  {
    if (type == (*iter).get_field_type()) 
    {
      MPI_Datatype dtype = (*iter).get_datatype();

      if (blocking_g)
      {
        send_recv((*iter).get_tx_ptr(), 
                  (*iter).get_rx_ptr(), 
                  dtype);
      }
      else
      {
        non_blocking((*iter), idx);
      }
    }

    ++iter;
  }

  num_used_ = idx;
}

void SubdomainBc::send_recv(const void *tx_ptr, void *rx_ptr, 
                            MPI_Datatype &t)
{
  MPI_Status status;
  int err;

  // Classic method
  if (rank_ > neighbour_) {
    err = MPI_Send(const_cast<void *>(tx_ptr), 1, t, neighbour_, 1, 
                   MPI_COMM_PHRED);
    err = MPI_Recv(rx_ptr, 1, t, neighbour_, 1, MPI_COMM_PHRED, &status);
  } else {
    err = MPI_Recv(rx_ptr, 1, t, neighbour_, 1, MPI_COMM_PHRED, &status);
    err = MPI_Send(const_cast<void *>(tx_ptr), 1, t, neighbour_, 1, 
                   MPI_COMM_PHRED);
  }
}

void SubdomainBc::non_blocking(RxTxData &data, int &idx)
{
  const void *tx_ptr = data.get_tx_ptr();
  void *rx_ptr = data.get_rx_ptr();
  MPI_Datatype t = data.get_datatype();
  int err;

  //req_[idx++] = data.get_send_req();
  //req_[idx++] = data.get_recv_req();

  if (rank_ > neighbour_) {
    err = MPI_Isend(const_cast<void *>(tx_ptr), 1, t, neighbour_, 1, 
                    MPI_COMM_PHRED, &req_[idx++]);

    if (err != MPI_SUCCESS)
      cerr << "ERROR ON RANK " << MPI_RANK << " calling MPI_Isend." << endl;

    err = MPI_Irecv(rx_ptr, 1, t, neighbour_, 1, MPI_COMM_PHRED, 
                    &req_[idx++]);

    if (err != MPI_SUCCESS)
      cerr << "ERROR ON RANK " << MPI_RANK << " calling MPI_Irecv." << endl;
  } else {
    err = MPI_Irecv(rx_ptr, 1, t, neighbour_, 1, MPI_COMM_PHRED, 
                    &req_[idx++]);

    if (err != MPI_SUCCESS)
      cerr << "ERROR ON RANK " << MPI_RANK << " calling MPI_Irecv." << endl;

    err = MPI_Isend(const_cast<void *>(tx_ptr), 1, t, neighbour_, 1, 
                    MPI_COMM_PHRED, &req_[idx++]);

    if (err != MPI_SUCCESS)
      cerr << "ERROR ON RANK " << MPI_RANK << " calling MPI_Isend." << endl;
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

  if (!blocking_g)
  {
    if (req_)
      delete[] req_;

    num_ = 2 * rx_tx_data_.size();
    req_ = new MPI_Request[num_];

    //memset(req_, 0, sizeof(MPI_Request) * num_);
  }
}

BoundaryCondition SubdomainBc::get_type() const
{
  return SUBDOMAIN;
}

void SubdomainBc::wait()
{
  if (!blocking_g)
    MPI_Waitall(num_used_, req_, MPI_STATUSES_IGNORE);
}
