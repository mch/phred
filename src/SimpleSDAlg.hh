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

#ifndef SIMPLE_SD_ALG_H
#define SIMPLE_SD_ALG_H

#include "SubdomainAlg.hh"

#include <math.h>

#include <iostream>

using namespace std;

/**
 * Implements a simple domain decomposition algorithm. 
 */
class SimpleSDAlg : public SubdomainAlg
{
private:
protected:
public:
  SimpleSDAlg();
  virtual ~SimpleSDAlg();

  /**
   * Implements a simple domain decomposition algorithm which divides
   * the domain into an even number (or possibly 3) of blocks, one of
   * which belongs to each processor. 
   *
   * @param rank the rank of the processor we are finding a grid for. 
   * @param size the total number of processors available to us. 
   * @param grid_info an object containing information about the
   * global grid as determined by parsing the input file. 
   *
   * @return a Grid object (class instance) which has its sizes set
   * but which has not yet allocated any memory.
   */
  GridInfo decompose_domain(GridInfo &info);  

};

#endif // SIMPLE_SD_ALG_H
