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

  /**
   * Implements an algorithm for domain decomposition using MPI
   * functions. Also creats a MPI_Communicator to help map processes
   * to hardware topology.
   *
   * @param grid_info an object containing information about the
   * global grid as determined by parsing the input file. 
   *
   * @return a Grid object (class instance) which has its sizes set
   * but which has not yet allocated any memory.
   */
  GridInfo decompose_domain(GridInfo &info);

private:

};

#endif // MPI_SUBDOMAIN_ALG_H
