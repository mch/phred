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

#include "MPISubdomainAlg.hh"
#include "Globals.hh"

#include <iostream>
using namespace std;

MPISubdomainAlg::MPISubdomainAlg()
{}

MPISubdomainAlg::~MPISubdomainAlg()
{}

GridInfo MPISubdomainAlg::decompose_domain(GridInfo &info)
{
  // Set up subdomains
  int dims[3], periods[3];

  for (int i = 0; i < 3; i++)
  {
    dims[i] = 0;
    periods[i] = 0;
  }

  MPI_Dims_create(MPI_SIZE, 3, dims);

  // Re-arrange the dims so that the largest number of MPI nodes are
  // placed along the longest grid axis.
  unsigned int sdx, sdy, sdz;
  sdx = info.global_dimx_;
  sdy = info.global_dimy_;
  sdz = info.global_dimz_;

  // sucktacular:
//   if (sdy > sdx)
//   {
//     int temp = dims[0];
//     dims[0] = dims[1];
//     dims[1] = temp;

//     if (sdz > sdy)
//     {
//       temp = dims[1];
//       dims[1] = dims[2];
//       dims[2] = temp;
//     }
//   } 
//   else if (sdz > sdx)
//   {
//     int temp = dims[0];
//     dims[0] = dims[2];
//     dims[2] = temp;
    
//     if (sdy > sdz)
//     {
//       temp = dims[1];
//       dims[1] = dims[2];
//       dims[2] = temp;
//     }
//   }

  // Set up a new communicator
  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &MPI_COMM_PHRED);

  // Create a new GridInfo for this rank, depending on where it is
  // within the topology.
  int coords[3], new_rank;
  MPI_Comm_rank(MPI_COMM_PHRED, &new_rank);
  MPI_Cart_coords(MPI_COMM_PHRED, new_rank, 3, coords);

#ifdef DEBUG
  cerr << "MPISubdomainAlg: old rank is " << MPI_RANK << ", new rank is " 
       << new_rank << endl;

  cerr << "This rank is at (" << coords[0] << ", " << coords[1] 
       << ", " << coords[2] << ") in the cartesian topology." << endl;
#endif

  MPI_RANK = new_rank;

  GridInfo result = info; 
  int n, m, p;

  n = dims[0];
  m = dims[1];
  p = dims[2];

  // Assign sizes and starting points including overlap
  result.dimx_ = static_cast<unsigned int>
    (floor(static_cast<double>(result.global_dimx_) / n));
  result.dimy_ = static_cast<unsigned int>
    (floor(static_cast<double>(result.global_dimy_) / m));
  result.dimz_ = static_cast<unsigned int>
    (floor(static_cast<double>(result.global_dimz_) / p));

  // In case the number of cells isn't evenly divisible
  int x = result.global_dimx_ % n;
  int y = result.global_dimy_ % m;
  int z = result.global_dimz_ % p;

  if (coords[0] == 0 && x != 0) {
    result.dimx_ += x;
  }

  if (coords[1] == 0 && y != 0) {
    result.dimy_ += y;
  }

  if (coords[2] == 0 && z != 0) {
    result.dimz_ += z;
  }

  result.start_x_no_sd_ = result.start_x_ = x + coords[0] * result.dimx_;
  result.start_y_no_sd_ = result.start_y_ = y + coords[1] * result.dimy_;
  result.start_z_no_sd_ = result.start_z_ = z * coords[2] * result.dimz_;
  result.dimx_no_sd_ = result.dimx_;
  result.dimy_no_sd_ = result.dimy_;
  result.dimz_no_sd_ = result.dimz_;

  SubdomainBc *sdbc = 0;

  int n_coords[3];
  int n_rank;

  if (coords[0] != 0) { // Back
    sdbc = new SubdomainBc();

    n_coords[0] = coords[0] - 1;
    n_coords[1] = coords[1];
    n_coords[2] = coords[2];
    MPI_Cart_rank(MPI_COMM_PHRED, n_coords, &n_rank);

    sdbc->set_neighbour(n_rank);
    sdbc->set_rank(MPI_RANK);
    result.set_boundary(BACK, sdbc);

    result.dimx_++;
    result.start_x_--;
  }

  if (coords[0] != dims[0] - 1) { // front
    sdbc = new SubdomainBc();

    n_coords[0] = coords[0] + 1;
    n_coords[1] = coords[1];
    n_coords[2] = coords[2];
    MPI_Cart_rank(MPI_COMM_PHRED, n_coords, &n_rank);

    sdbc->set_neighbour(n_rank);
    sdbc->set_rank(MPI_RANK);
    result.set_boundary(FRONT, sdbc);

    result.dimx_++;
  }

  if (coords[1] != 0) { // LEFT
    sdbc = new SubdomainBc();

    n_coords[0] = coords[0];
    n_coords[1] = coords[1] - 1;
    n_coords[2] = coords[2];
    MPI_Cart_rank(MPI_COMM_PHRED, n_coords, &n_rank);

    sdbc->set_neighbour(n_rank);
    sdbc->set_rank(MPI_RANK);
    result.set_boundary(LEFT, sdbc);

    result.dimy_++;
    result.start_y_--;
  }

  if (coords[1] != dims[1] - 1) { // RIGHT
    sdbc = new SubdomainBc();

    n_coords[0] = coords[0];
    n_coords[1] = coords[1] + 1;
    n_coords[2] = coords[2];
    MPI_Cart_rank(MPI_COMM_PHRED, n_coords, &n_rank);

    sdbc->set_neighbour(n_rank);
    sdbc->set_rank(MPI_RANK);
    result.set_boundary(RIGHT, sdbc);

    result.dimy_++;
  }

  if (coords[2] != 0) { // BOTTOM
    sdbc = new SubdomainBc();

    n_coords[0] = coords[0];
    n_coords[1] = coords[1];
    n_coords[2] = coords[2] - 1;
    MPI_Cart_rank(MPI_COMM_PHRED, n_coords, &n_rank);

    sdbc->set_neighbour(n_rank);
    sdbc->set_rank(MPI_RANK);
    result.set_boundary(BOTTOM, sdbc);

    result.dimz_++;
    result.start_z_--;
  }

  if (coords[2] != dims[2] - 1) { // TOP
    sdbc = new SubdomainBc();

    n_coords[0] = coords[0];
    n_coords[1] = coords[1];
    n_coords[2] = coords[2] + 1;
    MPI_Cart_rank(MPI_COMM_PHRED, n_coords, &n_rank);

    sdbc->set_neighbour(n_rank);
    sdbc->set_rank(MPI_RANK);
    result.set_boundary(TOP, sdbc);

    result.dimz_++;
  }

  return result;
}
