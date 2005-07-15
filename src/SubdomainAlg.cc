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

#include "SubdomainAlg.hh"
#include "Globals.hh"

#include <iostream>

using namespace std;

SubdomainAlg::SubdomainAlg()
{}

SubdomainAlg::~SubdomainAlg()
{}

GridInfo SubdomainAlg::decompose_domain(GridInfo &info)
{
  int dims[3], periods[3];;

  for (int i = 0; i < 3; i++)
  {
    dims[i] = 0;
    periods[i] = 0;
  }

  assign_processes(info, MPI_SIZE, dims);
  
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

  GridInfo result = calc_subdomain(info, coords, dims);

  setup_boundaries(result, coords, dims);

  return result;
}

void SubdomainAlg::setup_boundaries(GridInfo &gi, int coords[3], 
                                    int dims[3])
{
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
    gi.set_boundary(BACK, sdbc);

    gi.dimx_++;
    gi.start_x_--;

#ifdef DEBUG
    cerr << "MPISubdomainAlg: Rank " << MPI_RANK << " sharing data with "
         << n_rank << " at it's back face. " << endl;
#endif
  }

  if (coords[0] != dims[0] - 1) { // front
    sdbc = new SubdomainBc();

    n_coords[0] = coords[0] + 1;
    n_coords[1] = coords[1];
    n_coords[2] = coords[2];
    MPI_Cart_rank(MPI_COMM_PHRED, n_coords, &n_rank);

    sdbc->set_neighbour(n_rank);
    sdbc->set_rank(MPI_RANK);
    gi.set_boundary(FRONT, sdbc);

    gi.dimx_++;

#ifdef DEBUG
    cerr << "MPISubdomainAlg: Rank " << MPI_RANK << " sharing data with "
         << n_rank << " at it's front face. " << endl;
#endif
  }

  if (coords[1] != 0) { // LEFT
    sdbc = new SubdomainBc();

    n_coords[0] = coords[0];
    n_coords[1] = coords[1] - 1;
    n_coords[2] = coords[2];
    MPI_Cart_rank(MPI_COMM_PHRED, n_coords, &n_rank);

    sdbc->set_neighbour(n_rank);
    sdbc->set_rank(MPI_RANK);
    gi.set_boundary(LEFT, sdbc);

    gi.dimy_++;
    gi.start_y_--;

#ifdef DEBUG
    cerr << "MPISubdomainAlg: Rank " << MPI_RANK << " sharing data with "
         << n_rank << " at it's left face. " << endl;
#endif
  }

  if (coords[1] != dims[1] - 1) { // RIGHT
    sdbc = new SubdomainBc();

    n_coords[0] = coords[0];
    n_coords[1] = coords[1] + 1;
    n_coords[2] = coords[2];
    MPI_Cart_rank(MPI_COMM_PHRED, n_coords, &n_rank);

    sdbc->set_neighbour(n_rank);
    sdbc->set_rank(MPI_RANK);
    gi.set_boundary(RIGHT, sdbc);

    gi.dimy_++;

#ifdef DEBUG
    cerr << "MPISubdomainAlg: Rank " << MPI_RANK << " sharing data with "
         << n_rank << " at it's right face. " << endl;
#endif
  }

  if (coords[2] != 0) { // BOTTOM
    sdbc = new SubdomainBc();

    n_coords[0] = coords[0];
    n_coords[1] = coords[1];
    n_coords[2] = coords[2] - 1;
    MPI_Cart_rank(MPI_COMM_PHRED, n_coords, &n_rank);

    sdbc->set_neighbour(n_rank);
    sdbc->set_rank(MPI_RANK);
    gi.set_boundary(BOTTOM, sdbc);

    gi.dimz_++;
    gi.start_z_--;

#ifdef DEBUG
    cerr << "MPISubdomainAlg: Rank " << MPI_RANK << " sharing data with "
         << n_rank << " at it's bottom face. " << endl;
#endif
  }

  if (coords[2] != dims[2] - 1) { // TOP
    sdbc = new SubdomainBc();

    n_coords[0] = coords[0];
    n_coords[1] = coords[1];
    n_coords[2] = coords[2] + 1;
    MPI_Cart_rank(MPI_COMM_PHRED, n_coords, &n_rank);

    sdbc->set_neighbour(n_rank);
    sdbc->set_rank(MPI_RANK);
    gi.set_boundary(TOP, sdbc);

    gi.dimz_++;

#ifdef DEBUG
    cerr << "MPISubdomainAlg: Rank " << MPI_RANK << " sharing data with "
         << n_rank << " at it's top face. " << endl;
#endif
  }
}
