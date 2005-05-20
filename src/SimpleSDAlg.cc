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

/**********************************************************************
 * WARNING: THIS FILE HAS BEEN DEPRECATED
 **********************************************************************/

#include "SimpleSDAlg.hh"
#include "Globals.hh"
#include <exception>

SimpleSDAlg::SimpleSDAlg()
{}

SimpleSDAlg::~SimpleSDAlg()
{}

GridInfo SimpleSDAlg::decompose_domain(GridInfo &info)
{
  bool divided = false;

  if (MPI_SIZE < 0) // that's just wrong
    throw std::exception();

  unsigned int sz = static_cast<unsigned int>(MPI_SIZE);
  
  // n, m, and p are the number of divisions along the x, y, and z
  // axis' respectivly
  unsigned int sdx, sdy, sdz, n, m, p;
  sdx = info.global_dimx_;
  sdy = info.global_dimy_;
  sdz = info.global_dimz_;

  n = m = p = 1;

  for (int i = 0; i < MPI_SIZE; i++)
  {
    if (sdx >= sdy && sdx >= sdz && (n+1)*m*p <= sz) {
      n++;
      sdx = info.global_dimx_ / n;
      divided = true;
    }

    else if (sdy >= sdx && sdy >= sdz && n*(m+1)*p <= sz) {
      m++;
      sdy = info.global_dimy_ / m;
      divided = true;
    }

    else if (sdz >= sdx && sdz >= sdy && n*m*(p+1) <= sz) {
      p++;
      sdz = info.global_dimz_ / p;
      divided = true;
    }

    if (!divided) {
      if ((n+1)*m*p <= sz) {
        n++;
        sdx = info.global_dimx_ / n;
      }
      else if (n*(m+1)*p <= sz) {
        m++;
        sdy = info.global_dimy_ / m;
      }
      else if (n*m*(p+1) <= sz) {
        p++;
        sdz = info.global_dimz_ / p;
      }
    }

    divided = false;
  }

  if (n*m*p != sz) {
    cerr << "WARNING: simple domain decomposition did not result in one grid per\nprocessor!" << endl;
  }

  // Result initially has the exact same properties as the input.
  GridInfo result = info; 

  // Find the domain this grid will be in. 
  unsigned int x, y, z; 
  x = MPI_RANK % n;
  y = ((MPI_RANK - x)/n) % m;
  z = ((MPI_RANK - x)/n - y) / m;

  // Assign sizes and starting points including overlap
  result.dimx_ = static_cast<unsigned int>
    (floor(static_cast<double>(result.global_dimx_) / n));
  result.dimy_ = static_cast<unsigned int>
    (floor(static_cast<double>(result.global_dimy_) / m));
  result.dimz_ = static_cast<unsigned int>
    (floor(static_cast<double>(result.global_dimz_) / p));

  if (x == n-1) { // in case the number of cells isn't evenly divisible
    result.dimx_ = result.global_dimx_ - (n - 1)
      * static_cast<unsigned int>
      (floor(static_cast<double>(result.global_dimx_) / n));
  } 
  
  if (y == m-1) {
    result.dimy_ = result.global_dimy_ - (m - 1)
      * static_cast<unsigned int>
      (floor(static_cast<double>(result.global_dimy_) / m));
  } 

  if (z == p-1) {
    result.dimz_ = result.global_dimz_ - (p - 1)
      * static_cast<unsigned int>
      (floor(static_cast<double>(result.global_dimz_) / p));
  } 

  result.start_x_no_sd_ = result.start_x_ = x * result.dimx_;
  result.start_y_no_sd_ = result.start_y_ = y * result.dimy_;
  result.start_z_no_sd_ = result.start_z_ = z * result.dimz_;
  result.dimx_no_sd_ = result.dimx_;
  result.dimy_no_sd_ = result.dimy_;
  result.dimz_no_sd_ = result.dimz_;

  // Assign boundary conditions the ranks to talk to 
  SubdomainBc *sdbc = 0;

  if (x != 0) { // BACK
    sdbc = new SubdomainBc();

    sdbc->set_neighbour((z*m + y) * n + (x-1));
    sdbc->set_rank(MPI_RANK);
    result.set_boundary(BACK, sdbc);

    result.dimx_++;
    result.start_x_--;
  }

  if (x != n - 1) { // FRONT
    sdbc = new SubdomainBc();

    sdbc->set_neighbour((z*m + y) * n + (x+1));
    sdbc->set_rank(MPI_RANK);
    result.set_boundary(FRONT, sdbc);
    result.dimx_++;
  } 

  if (y != 0) { // LEFT 
    sdbc = new SubdomainBc();

    sdbc->set_neighbour((z*m + (y-1)) * n + x);
    sdbc->set_rank(MPI_RANK);
    result.set_boundary(LEFT, sdbc);
    result.dimy_++;
    result.start_y_--;
  } 

  if (y != m - 1) { // RIGHT
    sdbc = new SubdomainBc();

    sdbc->set_neighbour((z*m + (y+1)) * n + x);
    sdbc->set_rank(MPI_RANK);
    result.set_boundary(RIGHT, sdbc);
    result.dimy_++;
  }

  if (z != 0) { // BOTTOM
    sdbc = new SubdomainBc();

    sdbc->set_neighbour(((z-1)*m + y) * n + x);
    sdbc->set_rank(MPI_RANK);
    result.set_boundary(BOTTOM, sdbc);
    result.dimz_++;
    result.start_z_--;
  }

  if (z != p - 1) { // TOP
    sdbc = new SubdomainBc();

    sdbc->set_neighbour(((z+1)*m + y) * n + x);
    sdbc->set_rank(MPI_RANK);
    result.set_boundary(TOP, sdbc);
    result.dimz_++;
  }

  return result;
}
