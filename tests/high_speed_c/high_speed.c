/**
 * This file is a simple high speed test. The idea is to find out if
 * there is any significant impact resulting from the use of the C++
 * compiler to build phred, vs. using just a C compiler to build this
 * program. Of course, no C++ features are used here, and this program
 * is not very versitile. 
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

#include "../../src/config.h"
#include "../../src/Types.hh"
#include "../../src/Constants.hh"

/****************************************************************
 * Global Variables!
 ****************************************************************/

/* Number of materials we know about (0 is PEC)*/
unsigned int num_materials_;

/* E Field Material Coefficients */
mat_coef_t *Ca_;
mat_coef_t *Cbx_;
mat_coef_t *Cby_;
mat_coef_t *Cbz_;

/* H Field Coefficients */
mat_coef_t *Da_;
mat_coef_t *Dbx_;
mat_coef_t *Dby_;
mat_coef_t *Dbz_;

field_t *ex_;
field_t *ey_;
field_t *ez_;
field_t *hx_;
field_t *hy_;
field_t *hz_;

/* The material for each point in the grid. This is an index into
 * the material arrays, Ca, Cbx, etc. */
unsigned int *material_;

/* Time and space steppings; the distance between each point in the
 * grid. */ 
delta_t deltax_;
delta_t deltay_;
delta_t deltaz_;
delta_t deltat_;

/* Size of the grid along each dimension */
unsigned int dimx_;
unsigned int dimy_;
unsigned int dimz_;

/****************************************************************
 * Function declarations
 ****************************************************************/

void update_ex();
void update_ey();
void update_ez();

void update_hx();
void update_hy();
void update_hz();

unsigned int pi(unsigned int x, unsigned int y, 
                unsigned int z);

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
    update_ex();
    update_ey();
    update_ez();

    update_hx();
    update_hy();
    update_hz();

    ey_[pi(100, 100, 100)] = gaussm(i, 500e12, 1, 300e12);
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

/* Straight out of Taflove. */
/* This is a slow version */
void update_ex() 
{
  unsigned int mid, idx, idx2, i, j, k;
  
  for (i = 0; i < dimx_; i++) {
    for (j = 1; j < dimy_; j++) {

      idx = pi(i, j, 1);
      idx2 = pi(i, j-1, 1);

      for (k = 1; k < dimz_; k++) {
        mid = material_[idx];

        ex_[idx] = Ca_[mid] * ex_[idx]
          + Cby_[mid] * (hz_[idx] - hz_[idx2])
          + Cbz_[mid] * (hy_[idx-1] - hy_[idx]);
        
        idx++;
        idx2++;
      }
    }
  }
}

/* Straight out of Taflove. */
/* Pointer arithmetic should be faster. */
void update_ey() 
{
  unsigned int mid, i, j, k, idx;
  field_t *ey, *hx, *hz1, *hz2;

  for (i = 1; i < dimx_; i++) {
    for (j = 0; j < dimy_; j++) {

      idx = pi(i, j, 1);
      ey = &(ey_[idx]);
      hx = &(hx_[idx]);
      hz1 = &(hz_[pi(i-1, j, 1)]);
      hz2 = &(hz_[idx]);

      for (k = 1; k < dimz_; k++) {
        mid = material_[idx];

        *ey = Ca_[mid] * *ey
          + Cbz_[mid] * (*hx - *(hx-1))
          + Cbx_[mid] * (*hz1 - *hz2);

        ey++;
        hx++;
        hz1++;
        hz2++;
        idx++;
      }
    }
  }

}

/* Straight out of Taflove. */
/* Pointer arithmetic plus OpenMP threading should be faster still! */
void update_ez() 
{
  unsigned int mid, i, j, k, idx;
  field_t *ez, *hy1, *hy2, *hx1, *hx2;
  
  for (i = 1; i < dimx_; i++) {
    for (j = 1; j < dimy_; j++) {

      idx = pi(i, j, 0);
      ez = &(ez_[idx]);
      hy1 = &(hy_[idx]);
      hy2 = &(hy_[pi(i-1, j, 0)]);
      hx1 = &(hx_[pi(i, j-1, 0)]);
      hx2 = &(hx_[idx]);

#ifdef USE_OPENMP
#pragma opt parallel for \
        private(i, j, k, idx, idx2, mid, ez, hy1, hy2, hx1, hx2) \
        shared(Ca_, Cbx_, Cby_) \
        schedule(static)
#endif
      for (k = 0; k < dimz_; k++) {
        mid = material_[pi(i, j, k)];

        *ez = Ca_[mid] * *ez
          + Cbx_[mid] * (*hy1 - *hy2)
          + Cby_[mid] * (*hx1 - *hx2);

        ez++;
        hy1++; hy2++; hx1++; hx2++;
        idx++;
      }
    }
  }
}

void update_hx()
{
  unsigned int mid, i, j, k, idx;
  field_t *hx, *ez1, *ez2, *ey;

  for (i = 0; i < dimx_; i++) {
    for (j = 0; j < dimy_ - 1; j++) {

      idx = pi(i, j, 0);
      hx = &(hx_[idx]);
      ez1 = &(ez_[idx]);
      ez2 = &(ez_[pi(i, j+1, 0)]);
      ey = &(ey_[idx]);

      for (k = 0; k < dimz_ - 1; k++) {
        mid = material_[idx];

        *hx = Da_[mid] * *hx
          + Dby_[mid] * (*ez1 - *ez2)
          + Dbz_[mid] * (*(ey+1) - *ey);

        hx++; idx++;
        ez1++; ez2++; ey++;
      }
    }
  }
}

void update_hy()
{
  unsigned int mid, i, j, k, idx;
  field_t *hy, *ex, *ez1, *ez2;

  for (i = 0; i < dimx_ - 1; i++) {
    for (j = 0; j < dimy_; j++) {

      idx = pi(i, j, 0);
      hy = &(hy_[idx]);
      ex = &(ex_[idx]);
      ez1 = &(ez_[pi(i+1, j, 0)]);
      ez2 = &(ez_[idx]);

      for (k = 0; k < dimz_ - 1; k++) {
        mid = material_[idx];

        *hy = Da_[mid] * *hy
          + Dbz_[mid] * (*ex - *(ex + 1))
          + Dbx_[mid] * (*ez1 - *ez2);        
        
        hy++; idx++;
        ex++; ez1++; ez2++;
      }
    }
  }
}

void update_hz()
{
  unsigned int mid, i, j, k, idx;
  field_t *hz1, *ey1, *ey2, *ex1, *ex2;

  for (i = 0; i < dimx_ - 1; i++) {
    for (j = 0; j < dimy_ - 1; j++) {

      idx = pi(i, j, 0);
      hz1 = &(hz_[idx]);
      ey1 = &(ey_[idx]);
      ey2 = &(ey_[pi(i+1, j, 0)]);
      ex1 = &(ex_[pi(i, j+1, 0)]);
      ex2 = &(ex_[idx]);

      for (k = 0; k < dimz_; k++) {
        mid = material_[idx];

        *hz1 = Da_[mid] * *hz1
          + Dbx_[mid] * (*ey1 - *ey2)
          + Dby_[mid] * (*ex1 - *ex2);

        hz1++; idx++;
        ey1++; ey2++;
        ex1++; ex2++;
      }
    }
  }
}

unsigned int pi(unsigned int x, unsigned int y, 
                unsigned int z)
{
  return z + (y + x*dimy_) * dimz_;
}
