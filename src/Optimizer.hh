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

#ifndef OPTIMIZER_H
#define OPTIMIZER_H

/**
 * This class takes defines a problem to be optimized. It must be able
 * to set up the grid based on it's paramters, it must tell the
 * optimizer how much to change each parameter in order to calculate
 * the gradient, etc. 
 *
 * Sub-class this to create a problem. 
 */ 
class OptimProblem {
public:

protected:

private:
};

/**
 * This class takes a problem to be optimized and tries to optimize
 * it. From a given starting point, the gradient is computed. A
 * quadratic fit is done to solve for the step size, since a line
 * search would take too long. 
 */ 
class Optimizer {
public:

protected:

private:
};

#endif //  OPTIMIZER_H
