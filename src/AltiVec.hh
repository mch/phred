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

#ifndef ALTIVEC_H
#define ALTIVEC_H

#include "config.h"

#ifdef USE_ALTIVEC

/** 
 * Helpers for using the PowerPC processor AltiVec engine
 */

#define ZERO ((vector float)(0))

/* Unions make it easy to access elements in the vector */

typedef union
{
  vector signed char vec;
  signed char elements[16];
} CharVector_t;

typedef union 
{
  vector unsigned short vec;
  unsigned short elements[8];
} UnsignedShortVector_t; 

typedef union 
{
  vector signed short vec;
  signed short elements[8];
} ShortVector_t; 

typedef union 
{
  vector signed int vec;
  signed int elements[4];
} IntVector_t; 

typedef union 
{
  vector unsigned int vec;
  unsigned int elements[4];
} UnsignedIntVector_t; 

typedef union 
{
  vector float vec;
  float elements[4];
} FloatVector_t; 

/* Vector unit doesn't do doubles. */

#endif /* USE_ALTIVEC */

#endif /* ALTIVEC_H */
