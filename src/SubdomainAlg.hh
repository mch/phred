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

#ifndef SUBDOMAIN_ALG_H
#define SUBDOMAIN_ALG_H

#include "GridInfo.hh"

/**
 * This is an abstract base class which defines an interface for any
 * algorithm that produces a grid for a particular node. All
 * processors are intended to run this algorithm. 
 *
 * This ABC and it's subclasses implement the Strategy pattern from pg
 * 315 of gof1995. The GridInfo object is the Context the individual
 * strategies work with.
 */
class SubdomainAlg
{
private:
protected:
public:
  SubdomainAlg() {}
  virtual ~SubdomainAlg() {}

  /**
   * Subclasses must override this method and implement an algorithm
   * for domain decomposition. 
   *
   * @param grid_info an object containing information about the
   * global grid as determined by parsing the input file. 
   *
   * @return a Grid object (class instance) which has its sizes set
   * but which has not yet allocated any memory.
   */
  virtual GridInfo decompose_domain(GridInfo &info) = 0;

  /**
   * Subclasses must override this method and implement an algorithm
   * for domain decomposition. The input is a Region object describing
   * a set of voxels or Yee cells. The output is a new region which
   * describes a set of voxels local only to the currently running
   * process.
   *
   * @param region an object containing information about the
   * global grid as determined by parsing the input file. 
   *
   * @return a new OverlapRegion object describing the voxels that
   * belong only to this process, including any overlapping cells. 
   */
  //virtual OverlapRegion decompose_domain(Region &info) = 0;
};

#endif // SUBDOMAIN_ALG_H
