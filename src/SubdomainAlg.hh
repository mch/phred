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

#ifndef SUBDOMAIN_ALG_H
#define SUBDOMAIN_ALG_H

#include "GridInfo.hh"

/**
 * This is an abstract base class which defines an interface for any
 * algorithm that produces a grid for a particular node. All
 * processors are intended to run this algorithm. 
 *
 * All subclasses MUST: 
 * 1) Divide the global comuptational domain up into appropriatly
 *    sized chunks for each rank.
 * 2) Set up SubdomainBC boundary conditions where the local domains
 *    meet and information must be exchanged between MPI ranks.
 *
 * This ABC and it's subclasses implement the Strategy pattern from pg
 * 315 of gof1995. The GridInfo object is the Context the individual
 * strategies work with.
 */
class SubdomainAlg
{
public:
  SubdomainAlg();
  virtual ~SubdomainAlg();

  /**
   * This function returns a new GridInfo object for the local
   * domain. The input is the global GridInfo.
   * 
   * Subclasses may override this method to provide more sophisticated
   * behaviour than the default. The extra computational load caused
   * by PML boundary conditions is not considered by the default
   * implementation for instance.
   *
   * @param grid_info an object containing information about the
   * global grid as determined by parsing the input file. 
   *
   * @return a Grid object (class instance) which has its sizes set
   * but which has not yet allocated any memory.
   */
  virtual GridInfo decompose_domain(GridInfo &info);

protected:
  /**
   * This function divides a set of N processors into n processors
   * along the x axis, m processors along the y axis, and p processors
   * along the z axis. n*m*p == N.
   *
   * Subclasses must override this function. 
   */ 
  virtual void assign_processes(const GridInfo &info, int N, int dims[3]) = 0;

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
                                  const int dims[3]) = 0;

  /**
   * Sets up the Subdomain boundary conditions and adjusts the size of
   * the domain slightly to account for ghost cells.
   *
   * @param gi The GridInfo object for the local sub-domain. 
   * @param coords The (x,y,z) coordinates of the process within the
   * Cartesian topology.
   * @param dims The number of processors along each axis. 
   */ 
  void setup_boundaries(GridInfo &gi, int coords[3], int dims[3]);

};

#endif // SUBDOMAIN_ALG_H
