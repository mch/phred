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

/* Defines a bunch of useful constants that will be handy */

/* Use only C comments, or the AIX C compiler with complain */ 

/* Prefer static const over #define. */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "Types.hh"

// Pi
static const field_t PI = 3.14159265358979323846264338327;

// Permittivity of free space
static const field_t EPS_0 = 8.854185336732027e-12;

// Permeability of free space
static const field_t MU_0 = 1.256637061435917e-06;

// ?? 
static const field_t ETA_0 = 3.767303662405271e+02;

// Speed of light in free space
static const field_t C = 2.997925e8;

// Freespace wave impedenace, sqrt(MU_0 / EPS_0)
static const field_t ETA = 3.276730366241e2;

#endif

