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
#include <time.h>

#include "common.h"

/* #undef USE_OPENMP */

#ifdef USE_OPENMP
#include <omp.h>
#endif

void alloc_grid();
void free_grid();

void load_materials(unsigned int num_materials, field_t *eps,
                    field_t *mu);
void free_material();

field_t gaussm(unsigned int time_step, field_t deltaf, 
               field_t alpha, field_t f0);

void run(unsigned int num_time_steps);

/****************************************************************
 * MAIN
 ****************************************************************/
int main(int argc, char **argv)
{
  field_t eps[2], mu[2];
  unsigned int num_time_steps = 100, i, j, k; 
  FILE *fields;
  unsigned int x, y, z, numthreads;
  time_t start, now;
  clock_t cpu_start, cpu_now;

  x = 20;
  y = 20;
  z = 15;

  deltax_ = deltay_ = deltaz_ = 18.75e-9;
  deltat_ = 36e-18;

  dimx_ = 100; dimy_ = 100; dimz_ = 200;

  eps[0] = 1;
  eps[1] = 2.2;
  mu[0] = 1;
  mu[1] = 1;

  printf("High speed C FDTD test program\n");
  printf("Grid is %ix%ix%i, running for %i time steps.\n\n",
         dimx_, dimy_, dimz_, num_time_steps);

  load_materials(2, eps, mu);
  alloc_grid();

  for (i = 0; i < 40; i++)
    for (j = 0; j < 40; j++)
      for (k = 0; k < 40; k++)
	material_[pi(i, j, k)] = 1;

  /*fields = fopen("fields.txt", "w");*/

#ifdef USE_OPENMP

  printf("Number of threads in team: %i\n", omp_get_num_threads());
  printf("Maximum number of threads in team: %i\n", omp_get_max_threads());
  printf("Number of processors: %i\n", omp_get_num_procs());
  printf("Current thread number: %i\n", omp_get_thread_num());
  printf("Dynamic thread adjustment? %s\n",
	 (omp_get_dynamic() ? "yes" : "no"));
  printf("In parallel? %s\n",
	 (omp_in_parallel() ? "yes" : "no"));
  printf("Nested parallism? %s\n", 
	 (omp_get_nested() ? "yes" : "no"));

  /* Test the OpenMP */
  for (numthreads = 1; numthreads <= omp_get_max_threads(); numthreads++)
    {
      omp_set_num_threads(numthreads);
      start = time(NULL);
      cpu_start = clock();
      
      run(num_time_steps);
      
      now = time(NULL);
      cpu_now = clock();
      
      printf("%i of %i threads took %f wall clock seconds, and %f cpu seconds.\n", 
	     numthreads, omp_get_max_threads(),
	     (double)(now - start), (double)(cpu_now - cpu_start) / (double)CLOCKS_PER_SEC);
    }
#else
  run(num_time_steps);
#endif

  /*fclose(fields);*/

  /* Clean up */
  free_grid();

  return 0;
}

/****************************************************************
 * Function definitions
 ****************************************************************/
void run(unsigned int num_time_steps)
{
  unsigned int i;

  /* Run loop */
  for (i = 0; i < num_time_steps; i++)
  {
    printf("High speed C, time step %i, source: %g\n", 
           i, gaussm(i, 500e12, 1, 300e12));

#ifdef USE_OPENMP
    omp_e_update();
#else
    restricted_e_update();
    /* e_update(); */
#endif

    ey_[pi(20, 20, 20)] = ey_[pi(20, 20, 20)] + gaussm(i, 500e12, 1, 300e12);

#ifdef USE_OPENMP
    omp_h_update();
#else
    restricted_h_update();
    /* h_update(); */
#endif

/*      fprintf(fields, "%i %g %g %g %g %g %g %g\n",  */
/*  	    i, i * deltat_, ex_[pi(x, y, z)], */
/*  	    ey_[pi(x, y, z)], ez_[pi(x, y, z)],  */
/*  	    hx_[pi(x, y, z)], hy_[pi(x, y, z)], */
/*  	    hz_[pi(x, y, z)]); */
  }

}

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
    
    unsigned int zbytes = dimz_ * sizeof(mat_coef_t);
    Ca_temp_ = malloc(zbytes);
    Cbx_temp_ = malloc(zbytes);
    Cby_temp_ = malloc(zbytes);
    Cbz_temp_ = malloc(zbytes);

    Da_temp_ = malloc(zbytes);
    Dbx_temp_ = malloc(zbytes);
    Dby_temp_ = malloc(zbytes);
    Dbz_temp_ = malloc(zbytes);

    material_ = malloc(sz);

    if (ex_ && ey_ && ez_ && hx_ && hy_ && hz_ && material_
        && Ca_temp_ && Cbx_temp_ && Cby_temp_ && Cbz_temp_
        && Da_temp_ && Dbx_temp_ && Dby_temp_ && Dbz_temp_)
    {
      memset(ex_, 0, sz);
      memset(ey_, 0, sz);
      memset(ez_, 0, sz);
      
      memset(hx_, 0, sz);
      memset(hy_, 0, sz);
      memset(hz_, 0, sz);

      memset(material_, 0, sz);

      memset(Ca_temp_, 0, zbytes);
      memset(Cbx_temp_, 0, zbytes);
      memset(Cby_temp_, 0, zbytes);
      memset(Cbz_temp_, 0, zbytes);

      memset(Da_temp_, 0, zbytes);
      memset(Dbx_temp_, 0, zbytes);
      memset(Dby_temp_, 0, zbytes);
      memset(Dbz_temp_, 0, zbytes);

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

    free(Ca_temp_);
    free(Cbx_temp_);
    free(Cby_temp_);
    free(Cbz_temp_);

    free(Da_temp_);
    free(Dbx_temp_);
    free(Dby_temp_);
    free(Dbz_temp_);
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



