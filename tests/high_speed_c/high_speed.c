/**
 * This file is a simple high speed test. The idea is to find out if
 * there is any significant impact resulting from the use of the C++
 * compiler to build phred, vs. using just a C compiler to build this
 * program. Of course, no C++ features are used here, and this program
 * is not very versatile. 
 *
 * This program is single node only, no MPI, maybe OpenMP and Altivec,
 * if it's worth the trouble to test those. 
 *
 * This is short on error handling, so don't mess it up. 
 */

/* Compile with "cc_r -pg -O2 -qSMP=omp -o hsfdtd -lm_r high_speed.c" on AIX */
/* Set OMP_NUM_THREADS to the number of processors for it to work. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "common.h"

void alloc_grid();
void free_grid();

void load_materials(unsigned int num_materials, field_t *eps,
                    field_t *mu);
void free_material();

field_t gaussm(unsigned int time_step, field_t deltaf, 
               field_t alpha, field_t f0);

/****************************************************************
 * MAIN
 ****************************************************************/
int main(int argc, char **argv)
{
  field_t eps[2], mu[2];
  unsigned int num_time_steps = 50, i; 

  deltax_ = deltay_ = deltaz_ = 18.75e-9;
  deltat_ = 36e-18;

  dimx_ = dimy_ = dimz_ = 250;

  eps[0] = 1;
  eps[1] = 2.2;
  mu[0] = 1;
  mu[1] = 1;

  printf("High speed C FDTD test program\n");
  printf("Grid is %ix%ix%i, running for %i time steps.\n\n",
         dimx_, dimy_, dimz_, num_time_steps);

  load_materials(2, eps, mu);
  alloc_grid();

  /* Run loop */
  for (i = 0; i < num_time_steps; i++)
  {
    printf("High speed C, time step %i, source: %g\n", 
           i, gaussm(i, 500e12, 1, 300e12));
    e_update();

    ey_[pi(100, 100, 100)] = gaussm(i, 500e12, 1, 300e12);

    h_update();

  }

  /* Clean up */
  free_grid();

  return 0;
}

/****************************************************************
 * Function definitions
 ****************************************************************/

field_t gaussm(unsigned int time_step, field_t deltaf, 
               field_t alpha, field_t f0)
{
  field_t t = (field_t)time_step * deltat_;

  return alpha * exp(-pow((t - 4. / (PI * deltaf)) * deltaf * PI, 2))
    * sin(2. * PI * f0 * (t - 4. / (PI * deltaf)));
}

void alloc_grid()
{
  unsigned int sz = 0;

  sz = dimx_ * dimy_ * dimz_ * sizeof(field_t);

  if (sz > 0) 
  {
    /* Fortunatly for us, on OS X malloc returns memory block aligned
       on 16 byte boundaries required by AltiVec */ 
    ex_ = malloc(sz);
    ey_ = malloc(sz);
    ez_ = malloc(sz);
    
    hx_ = malloc(sz);
    hy_ = malloc(sz);
    hz_ = malloc(sz);
    
    material_ = malloc(sz);

    if (ex_ && ey_ && ez_ && hx_ && hy_ && hz_ && material_)
    {
      memset(ex_, 0, sz);
      memset(ey_, 0, sz);
      memset(ez_, 0, sz);
      
      memset(hx_, 0, sz);
      memset(hy_, 0, sz);
      memset(hz_, 0, sz);

      memset(material_, 0, sz);
    } else {
      printf("Insufficient memory; try reducing the grid size. \n");
      exit(1);
    }
  }
}

void free_grid()
{
  if (ex_ || ey_ || ez_ || hx_ || hy_ || hz_) {
    free(ex_);
    free(ey_);
    free(ez_);

    free(hx_);
    free(hy_);
    free(hz_);

    free(material_);
  }
}

void load_materials(unsigned int num_materials, field_t *epsarg,
                    field_t *muarg)
{
  size_t sz = (num_materials + 1) * sizeof(mat_coef_t);
  unsigned int index = 0;
  mat_prop_t eps;
  mat_prop_t sig;
  mat_prop_t mu;
  mat_prop_t sigs;

  Ca_ = malloc(sz);
  Da_ = malloc(sz);

  Cbx_ = malloc(sz);
  Dbx_ = malloc(sz);

  /* Save some memory if possible. */
  if (deltay_ == deltax_)
  {
    Cby_ = Cbx_;
    Dby_ = Dbx_;
  } else {
    Cby_ = malloc(sz);
    Dby_ = malloc(sz);
  }

  if (deltaz_ == deltax_)
  {
    Cbz_ = Cbx_;
    Dbz_ = Dbx_;
  } else if (deltaz_ == deltay_) {
    Cbz_ = Cby_;
    Dbz_ = Dby_;
  } else {
    Cbz_ = malloc(sz);
    Dbz_ = malloc(sz);
  }

  /* The first one is always PEC */
  Ca_[index] = 1;
  Cbx_[index] = Cby_[index] = Cbz_[index] = 0;

  Da_[index] = 1;
  Dbx_[index] = Dby_[index] = Dbz_[index] = 0;
  
  for (index = 1; index < num_materials; index++)
  {
    eps = epsarg[index - 1] * EPS_0;
    sig = 0;
    mu = muarg[index - 1] * MU_0;
    sigs = 0;

    Ca_[index] = (1 - (sig * deltat_ * 0.5)/eps) / 
                 (1 + (sig * deltat_ * 0.5)/eps);

    Da_[index] = (1 - (sigs * deltat_ * 0.5)/mu) / 
                 (1 + (sigs * deltat_ * 0.5)/mu);

    
    Cbx_[index] = (deltat_ / (eps * deltax_)) / 
                  (1 + (sig * deltat_ * 0.5)/eps);

    Dbx_[index] = (deltat_ / (mu * deltax_)) / 
                  (1 + (sigs * deltat_ * 0.5)/mu);

    if (deltay_ != deltax_)
    {    
      Cby_[index] = (deltat_ / (eps * deltay_)) / 
                    (1 + (sig * deltat_ * 0.5)/eps);

      Dby_[index] = (deltat_ / (mu * deltay_)) / 
                    (1 + (sigs * deltat_ * 0.5)/mu);
    }

    if (deltaz_ != deltax_ && deltaz_ != deltay_)
    {
      Cbz_[index] = (deltat_ / (eps * deltaz_)) / 
                    (1 + (sig * deltat_ * 0.5)/eps);

      Dbz_[index] = (deltat_ / (mu * deltaz_)) / 
                    (1 + (sigs * deltat_ * 0.5)/mu);
    }
  }
}

void free_material()
{
  free(Ca_);
  free(Da_);

  free(Cbx_);
  free(Dbx_);

  /* Save some memory if possible. */
  if (deltay_ != deltax_)
  {
    free(Cby_);
    free(Dby_);
  }

  if (deltaz_ != deltax_ && deltaz_ != deltay_) {
    free(Cbz_);
    free(Dbz_);
  }
}

unsigned int pi(unsigned int x, unsigned int y, 
                unsigned int z)
{
  return z + (y + x*dimy_) * dimz_;
}



