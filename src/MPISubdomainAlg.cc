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

#include <cmath>
#include <iostream>
using namespace std;

MPISubdomainAlg::MPISubdomainAlg()
{}

MPISubdomainAlg::~MPISubdomainAlg()
{}

void MPISubdomainAlg::assign_processes(const GridInfo &info, int sz, 
                                       int dims[3])
{
  unsigned int sdx, sdy, sdz, n, m, p;
  bool divided = false;
  
  // The following divides up the domain such that the message size is
  // minimized. n, m, and p are the number of divisions along the x,
  // y, and z axis' respectivly

  sdx = info.global_dimx_;
  sdy = info.global_dimy_;
  sdz = info.global_dimz_;
    
  n = m = p = 1;
    
  for (int i = 0; i < MPI_SIZE; i++)
  {

    // Add a weighting factor here to express the fact that data
    // transmission is faster when contiguous memory is involved, such
    // as when the YZ plane must be transmitted (domain divided along X
    // axis).
    
    if (sdx >= sdy && /*sdx >= sdz &&*/ (n+1)*m*p <= sz) {
      n++;
      sdx = info.global_dimx_ / n;
      divided = true;
    }
      
    else if (sdy >= sdx /*&& sdy >= sdz*/ && n*(m+1)*p <= sz) {
      m++;
      sdy = info.global_dimy_ / m;
      divided = true;
    }
      
    // TEMPORARY, UNTIL UPML with z divisions is fixed
//     else if (sdz >= sdx && sdz >= sdy && n*m*(p+1) <= sz) {
//       p++;
//       sdz = info.global_dimz_ / p;
//       divided = true;
//     }
      
    if (!divided) {
      if ((n+1)*m*p <= sz) {
        n++;
        sdx = info.global_dimx_ / n;
      }
      else if (n*(m+1)*p <= sz) {
        m++;
        sdy = info.global_dimy_ / m;
      }

      // TEMPORARY, UNTIL UPML with z divisions is fixed
//       else if (n*m*(p+1) <= sz) {
//         p++;
//         sdz = info.global_dimz_ / p;
//       }
    }

    divided = false;
  }
    
  if (n*m*p != sz) {
#ifdef DEBUG
    cerr << "MPISubdomainAlg: Algorithm unable to perform decomposition. \n" 
         << "Using MPI instead. " << endl;
#endif
    // This uses MPI to divide up the domain, but it doesn't know
    // anything about minimizing message size, or other factors. The
    // good thing is that it handles any MPI_SIZE.
    
    // TEMPORARY, UNTIL UPML with z divisions is fixed
    dims[2] = 1;

    //  MPI_DIMS_CREATE chooses dimensions so that the resulting grid
    //  is as close as possible to being an ndims-dimensional cube.
    MPI_Dims_create(MPI_SIZE, 3, dims);
    
    n = dims[0];
    m = dims[1];
    p = dims[2];
  }
  
  dims[0] = n;
  dims[1] = m;
  dims[2] = p;
}

GridInfo MPISubdomainAlg::calc_subdomain(const GridInfo &gi,
                                         const int coords[3], 
                                         const int dims[3])
{
  GridInfo result = gi; 
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
  int x = result.global_dimx_ - n * result.dimx_;
  int y = result.global_dimy_ - m * result.dimy_;
  int z = result.global_dimz_ - p * result.dimz_;

  result.start_x_no_sd_ = result.start_x_ = x + coords[0] * result.dimx_;
  result.start_y_no_sd_ = result.start_y_ = y + coords[1] * result.dimy_;
  result.start_z_no_sd_ = result.start_z_ = z + coords[2] * result.dimz_;

  if (coords[0] < x) {
    result.dimx_ += 1;
    
    result.start_x_no_sd_ = result.start_x_ 
      = result.start_x_ - (x - coords[0]);
  }

  if (coords[1] < y) {
    result.dimy_ += 1;

    result.start_y_no_sd_ = result.start_y_ 
      = result.start_y_ - (y - coords[1]);
  }

  if (coords[2] < z) {
    result.dimz_ += 1;

    result.start_z_no_sd_ = result.start_z_ 
      = result.start_z_ - (z - coords[2]);
  }

  result.dimx_no_sd_ = result.dimx_;
  result.dimy_no_sd_ = result.dimy_;
  result.dimz_no_sd_ = result.dimz_;
  
  return result;
}
