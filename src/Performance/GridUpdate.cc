/* 
   Phred - Phred is a parallel finite difference time domain
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

#ifndef GRID_UPDATE_H
#define GRID_UPDATE_H

// All of the following updates use a form of flux density
// formulation, to allow for superposition of two seperate updates. 

/**
 * Calculates the electric field intensity from the electric flux
 * density 
 */
class IntensityUpdate
{

};

/**
 * Non-lossy dielectric update. Fast and simple. 
 */ 
class DielectricUpdate
{
  
};

/**
 * Lossy dielectric update.
 */ 
class LossyDielectricUpdate
{

};

/**
 * Unmagnetized plasma update.
 */
class PlasmaUpdate
{

};

/**
 * PML Update
 */ 
class PmlUpdate
{

};


/**
 * This class contains templatable loops for updating a region of the
 * grid. The idea to to pull as much of the looping constructs and set
 * up code as possible back into this module, so that modules which
 * implement the actual updates and dispersions are boiled down to
 * their most simple essences.
 */ 
template<class T> // Use the Boost Preprocessor libs to make this
                  // class accept multiple args...
class GridUpdate
{
public:


private:

};

#endif // GRID_UPDATE_H
