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

#include "Periodic.hh"
#include "../Globals.hh" 

Periodic::Periodic(shared_ptr<PeriodicExcitation> pe)
  : pe_(pe), valid_(false)
{
  for (int i = 0; i < 6; i++)
  {
    faces_[i] = false;
    exchange_rank_[i] = MPI_RANK;
  }
}

Periodic::~Periodic()
{

}

void Periodic::apply(Face face, Grid &grid, FieldType type)
{
  //if (!valid_)
  //{
    //throw BoundaryConditionException("Periodic: boundaries are not valid.");
  //}

  if (faces_[face])
  {
    if (type == E) 
      copy_e(face, grid);

    if (type == H) 
      copy_h(face, grid);

  }
}

void Periodic::init(const Grid &grid, Face face)
{
  // This is where we find out what planes we've been applied to. 
  faces_[face] = true;

  bool fbv = false;
  bool lrv = false;
  bool tbv = false; 

  // These tests have to rely on the global face settings, pre domain decomp. 
  if (faces_[FRONT] && faces_[BACK]) 
    fbv = true; 

  if (faces_[LEFT] && faces_[RIGHT])
    lrv = true;

  if (faces_[TOP] && faces_[BOTTOM])
    tbv = true; 

  Face ex_face = pe_->get_face();

  if ((ex_face == FRONT || ex_face == BACK) && fbv)
    valid_ = false;

  else if ((ex_face == LEFT || ex_face == RIGHT) && lrv)
    valid_ = false;

  else if ((ex_face == TOP || ex_face == BOTTOM) && tbv)
    valid_ = false;

  else if (fbv || tbv || lrv)
    valid_ = true; 

  // Find out where we are
  int coords[3], dims[3], periods[3];
  MPI_Cart_get(MPI_COMM_PHRED, 3, dims, periods, coords);

  // Find out the rank with which data must be exchanged
  int n_coords[3];
  int n_rank;

  switch (face)
  {
  case BACK:
    n_coords[0] = dims[0] - 1;
    n_coords[1] = coords[1];
    n_coords[2] = coords[2];
    break;

  case FRONT:
    n_coords[0] = 0;
    n_coords[1] = coords[1];
    n_coords[2] = coords[2];
    break;

  case LEFT:
    n_coords[0] = coords[0];
    n_coords[1] = dims[1] - 1;
    n_coords[2] = coords[2];
    break;

  case RIGHT:
    n_coords[0] = coords[0];
    n_coords[1] = 0;
    n_coords[2] = coords[2];
    break;

  case BOTTOM:
    n_coords[0] = coords[0];
    n_coords[1] = coords[1];
    n_coords[2] = dims[2] - 1;
    break;

  case TOP:
    n_coords[0] = coords[0];
    n_coords[1] = coords[1];
    n_coords[2] = 0;
    break;
  }

  MPI_Cart_rank(MPI_COMM_PHRED, n_coords, &n_rank);

  exchange_rank_[face] = n_rank;
}

void Periodic::deinit(const Grid &grid, Face face)
{

}

void Periodic::copy_e(Face face, Grid &grid)
{
  MPI_Status status;
  int idx = grid.pi(grid.get_ldx_sd() - 1, 0, 0);
  int mpi_err = 0;
  int sz = 0;
  MPI_Aint extnt = 0;
  
  MPI_Datatype t; 
  t = grid.get_plane_dt(face);

  switch (face)
  {
  case BACK:
    if (MPI_RANK != exchange_rank_[face])
    {
      MPI_Recv(grid.get_ey_ptr(grid.pi(0, 0, 0)), 1, t, 
               exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
      MPI_Recv(grid.get_ez_ptr(grid.pi(0, 0, 0)), 1, t, 
               exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
    } 
    else
    {
      MPI_Sendrecv(grid.get_ey_ptr(idx), 1, 
                   t, MPI_RANK, 1, 
                   grid.get_ey_ptr(0), 1, 
                   t, exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
      
      MPI_Sendrecv(grid.get_ez_ptr(idx), 1, 
                   t, MPI_RANK, 1, 
                   grid.get_ez_ptr(0), 1, 
                   t, exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
    }
    break;
    
  case FRONT:
    if (MPI_RANK != exchange_rank_[face])
    {
      MPI_Send(grid.get_ey_ptr(grid.pi(grid.get_ldx_sd() - 1, 0, 0)), 1, 
               t, exchange_rank_[face], 1, MPI_COMM_PHRED);
      
      MPI_Send(grid.get_ez_ptr(grid.pi(grid.get_ldx_sd() - 1, 0, 0)), 1, 
               t, exchange_rank_[face], 1, MPI_COMM_PHRED);
    }
    
    break;

  case LEFT:
    if (MPI_RANK != exchange_rank_[face])
    {
      MPI_Recv(grid.get_ex_ptr(grid.pi(0, 0, 0)), 1, t, 
               exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
      
      MPI_Recv(grid.get_ez_ptr(grid.pi(0, 0, 0)), 1, t, 
               exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
    }
    else
    {
      MPI_Sendrecv(grid.get_ex_ptr(grid.pi(0, grid.get_ldy_sd() - 1, 0)), 1, 
                   t, MPI_RANK, 1, grid.get_ex_ptr(grid.pi(0, 0, 0)), 1, 
                   t, exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
      
      MPI_Sendrecv(grid.get_ez_ptr(grid.pi(0, grid.get_ldy_sd() - 1, 0)), 1, 
                   t, MPI_RANK, 1, grid.get_ez_ptr(grid.pi(0, 0, 0)), 1, 
                   t, exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
    }
    
    break;

  case RIGHT:
    if (MPI_RANK != exchange_rank_[face])
    {
      MPI_Send(grid.get_ex_ptr(grid.pi(0, grid.get_ldy_sd() - 1, 0)), 1, 
               t, exchange_rank_[face], 1, MPI_COMM_PHRED);
      
      MPI_Send(grid.get_ez_ptr(grid.pi(0, grid.get_ldy_sd() - 1, 0)), 1, 
               t, exchange_rank_[face], 1, MPI_COMM_PHRED);
    }

    break;

  case BOTTOM:
    if (MPI_RANK != exchange_rank_[face])
    {
      MPI_Recv(grid.get_ex_ptr(grid.pi(0, 0, 0)), 1, t, 
               exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
      
      MPI_Recv(grid.get_ey_ptr(grid.pi(0, 0, 0)), 1, t, 
               exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
    }
    else
    {
      MPI_Sendrecv(grid.get_ex_ptr(grid.pi(0, 0, grid.get_ldz_sd() - 1)), 1, 
                   t, MPI_RANK, 1, 
                   grid.get_ex_ptr(grid.pi(0, 0, 0)), 1, 
                   t, exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
      
      MPI_Sendrecv(grid.get_ey_ptr(grid.pi(0, 0, grid.get_ldz_sd() - 1)), 1, 
                   t, MPI_RANK, 1, 
                   grid.get_ey_ptr(grid.pi(0, 0, 0)), 1, 
                   t, exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
    }

    break;

  case TOP:
    if (MPI_RANK != exchange_rank_[face])
    {
      MPI_Send(grid.get_ex_ptr(grid.pi(0, 0, grid.get_ldz_sd() - 1)), 1, 
               t, exchange_rank_[face], 1, MPI_COMM_PHRED);
      
      MPI_Send(grid.get_ey_ptr(grid.pi(0, 0, grid.get_ldz_sd() - 1)), 1, 
               t, exchange_rank_[face], 1, MPI_COMM_PHRED);
    }

    break;
  }
  
}


void Periodic::copy_h(Face face, Grid &grid)
{
  MPI_Status status;
  MPI_Datatype t; 
  t = grid.get_plane_dt(face);

  switch (face)
  {
  case FRONT:
    if (MPI_RANK != exchange_rank_[face])
    {
      MPI_Recv(grid.get_hy_ptr(grid.pi(grid.get_ldx_sd() - 1, 0, 0)), 1, 
               t, exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
      
      MPI_Recv(grid.get_hz_ptr(grid.pi(grid.get_ldx_sd() - 1, 0, 0)), 1, 
               t, exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
    }
    else
    {
      MPI_Sendrecv(grid.get_hy_ptr(grid.pi(0, 0, 0)), 1, 
                   t, MPI_RANK, 1, 
                   grid.get_hy_ptr(grid.pi(grid.get_ldx_sd() - 1, 0, 0)), 1, 
                   t, exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
      
      MPI_Sendrecv(grid.get_hz_ptr(grid.pi(0, 0, 0)), 1, 
                   t, MPI_RANK, 1, 
                   grid.get_hz_ptr(grid.pi(grid.get_ldx_sd() - 1, 0, 0)), 1, 
                   t, exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
    }

   break;

  case BACK:
    if (MPI_RANK != exchange_rank_[face])
    {
      MPI_Send(grid.get_hy_ptr(grid.pi(0, 0, 0)), 1, 
               t, exchange_rank_[face], 1, MPI_COMM_PHRED);
      
      MPI_Send(grid.get_hz_ptr(grid.pi(0, 0, 0)), 1, 
               t, exchange_rank_[face], 1, MPI_COMM_PHRED);
    }

    break;

  case RIGHT:
    if (MPI_RANK != exchange_rank_[face])
    {
      MPI_Recv(grid.get_hx_ptr(grid.pi(0, grid.get_ldy_sd() - 1, 0)), 1, 
               t, exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
      
      MPI_Recv(grid.get_hz_ptr(grid.pi(0, grid.get_ldy_sd() - 1, 0)), 1, 
               t, exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
    } 
    else
    {
      MPI_Sendrecv(grid.get_hx_ptr(grid.pi(0, 0, 0)), 1, 
                   t, MPI_RANK, 1, 
                   grid.get_hx_ptr(grid.pi(0, grid.get_ldy_sd() - 1, 0)), 1, 
                   t, exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
      
      MPI_Sendrecv(grid.get_hz_ptr(grid.pi(0, 0, 0)), 1, 
                   t, MPI_RANK, 1, 
                   grid.get_hz_ptr(grid.pi(0, grid.get_ldy_sd() - 1, 0)), 1, 
                   t, exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
    }

    break;

  case LEFT:
    if (MPI_RANK != exchange_rank_[face])
    {
      MPI_Send(grid.get_hx_ptr(grid.pi(0, 0, 0)), 1, 
               t, exchange_rank_[face], 1, MPI_COMM_PHRED);
      
      MPI_Send(grid.get_hz_ptr(grid.pi(0, 0, 0)), 1, 
               t, exchange_rank_[face], 1, MPI_COMM_PHRED);
    }

    break;

  case TOP:
    if (MPI_RANK != exchange_rank_[face])
    {
      MPI_Recv(grid.get_hx_ptr(grid.pi(0, 0, grid.get_ldz_sd() - 1)), 1, 
               t, exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
      
      MPI_Recv(grid.get_hy_ptr(grid.pi(0, 0, grid.get_ldz_sd() - 1)), 1, 
               t, exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
    }
    else
    {
      MPI_Sendrecv(grid.get_hx_ptr(grid.pi(0, 0, 0)), 1, 
                   t, MPI_RANK, 1, 
                   grid.get_hx_ptr(grid.pi(0, 0, grid.get_ldz_sd() - 1)), 1, 
                   t, exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
      
      MPI_Sendrecv(grid.get_hy_ptr(grid.pi(0, 0, 0)), 1, 
                   t, MPI_RANK, 1, 
                   grid.get_hy_ptr(grid.pi(0, 0, grid.get_ldz_sd() - 1)), 1, 
                   t, exchange_rank_[face], 1, MPI_COMM_PHRED, &status);
    }

    break;

  case BOTTOM:

    if (MPI_RANK != exchange_rank_[face])
    {
      MPI_Send(grid.get_hx_ptr(grid.pi(0, 0, 0)), 1, 
               t, exchange_rank_[face], 1, MPI_COMM_PHRED);
      
      MPI_Send(grid.get_hy_ptr(grid.pi(0, 0, 0)), 1, 
               t, exchange_rank_[face], 1, MPI_COMM_PHRED);
    }

    break;
  }

}

