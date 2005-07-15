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

#ifndef MPI_SUBDOMAIN_ALG_H
#define MPI_SUBDOMAIN_ALG_H

#include "SubdomainAlg.hh"

/**
 * This class implements what I should have done in the first place:
 * it uses the MPI_Cart_Create and MPI_Dims_create functions to
 * compute the subdomain distribution. Minor adjustments are made for
 * ghost cells. 
 *
 * Unlike SimpleSubdomainAlg, this does not attempt to optimize the
 * domain decomposition in any way.
 */
class MPISubdomainAlg : public SubdomainAlg {
public:
  MPISubdomainAlg();
  ~MPISubdomainAlg();

protected:
  /**
   * This function divides a set of N processors into n processors
   * along the x axis, m processors along the y axis, and p processors
   * along the z axis. n*m*p == N.
   *
   * Subclasses must override this function. 
   */ 
  virtual void assign_processes(const GridInfo &info, int N, int dims[3]);

  /**
   * Calculate the size and start of the subdomain on a specific rank.
   *
   * @param gi GridInfo for global computational domain
   * @param coords The (x,y,z) coordinates of the process within the
   * Cartesian topology.
   * @param dims The number of processors along the x, y, and z axis. 
   * @return a new GridInfo object based on gi which represents the
   * subdomain on process rank.
   */ 
  virtual GridInfo calc_subdomain(const GridInfo &gi,
                                  const int coords[3], 
                                  const int dims[3]);

};

#endif // MPI_SUBDOMAIN_ALG_H
