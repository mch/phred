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
