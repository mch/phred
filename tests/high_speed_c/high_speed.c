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
#include <unistd.h>

/* For pointer arithmetic */ 
#include <inttypes.h>

#include "common.h"

/* #undef USE_OPENMP */

#ifdef USE_OPENMP
#include <omp.h>
#endif

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

/* Temporary coefficient holders. Pulls all the non-contiguous memory
   access out of the main update loop... hopefully making it faster. */
mat_coef_t * restrict Ca_temp_;
mat_coef_t * restrict Cbx_temp_;
mat_coef_t * restrict Cby_temp_;
mat_coef_t * restrict Cbz_temp_;

mat_coef_t * restrict Da_temp_;
mat_coef_t * restrict Dbx_temp_;
mat_coef_t * restrict Dby_temp_;
mat_coef_t * restrict Dbz_temp_;

mat_coef_t * restrict Ca_temp_orig_;
mat_coef_t * restrict Cbx_temp_orig_;
mat_coef_t * restrict Cby_temp_orig_;
mat_coef_t * restrict Cbz_temp_orig_;

mat_coef_t * restrict Da_temp_orig_;
mat_coef_t * restrict Dbx_temp_orig_;
mat_coef_t * restrict Dby_temp_orig_;
mat_coef_t * restrict Dbz_temp_orig_;

/* The fields that get operated on. These are padded a bit so that
   better use of cache can be made. */ 
field_t *ex_;
field_t *ey_;
field_t *ez_;
field_t *hx_;
field_t *hy_;
field_t *hz_;

/* These are the origionally allocated pointers for field data. Used
   to free() the memory. */ 
field_t *ex_orig_;
field_t *ey_orig_;
field_t *ez_orig_;
field_t *hx_orig_;
field_t *hy_orig_;
field_t *hz_orig_;

/* Temporary data if needed */ 
field_t *temp;

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

/* Padding along the z and y axis such that the length of a block of
   memory along the z axis is a multiple of 16 bytes. */
int z_padding_;
int y_padding_;

/****************************************************************
 * Functions for this file
 ****************************************************************/
void alloc_grid();
void free_grid();

void load_materials(unsigned int num_materials, field_t *eps,
                    field_t *mu);
void free_material();

field_t gaussm(unsigned int time_step, field_t deltaf, 
               field_t alpha, field_t f0);

void run(unsigned int num_time_steps);

/****************************************************************
 * Da switches
 ****************************************************************/
char cache_padding_g_;
char quiet_g_;
char openmp_loop_g_;
char simple_loop_g_;
char cache_loop_g_;
char combined_loop_g_;
char combined_simple_mat_loop_g_;
int x_size_g_;
int y_size_g_;
int z_size_g_;
int n_g_;

static int
decode_switches (int argc, char **argv)
{
  int c;
  char *arg = 0;
  opterr = 0;

  while ((c = getopt (argc, argv, "mqoscCx:y:z:n:")) != -1)
    switch (c)
    {
    case 'o':
#ifdef USE_OPENMP
      openmp_loop_g_ = 1;
      simple_loop_g_ = 0;
      cache_loop_g_ = 0;
      combined_loop_g_ = 0;
      combined_simple_mat_loop_g_ = 0;
#else
      printf("OpenMP is not available.\n");
#endif
      break;

    case 'q':
      quiet_g_ = 1;
      break;

    case 's':
      simple_loop_g_ = 1;
      openmp_loop_g_ = 0;
      cache_loop_g_ = 0;
      combined_loop_g_ = 0;
      combined_simple_mat_loop_g_ = 0;
      break;

    case 'C':
      cache_loop_g_ = 1;
      simple_loop_g_ = 0;
      openmp_loop_g_ = 0;
      combined_loop_g_ = 0;
      combined_simple_mat_loop_g_ = 0;
      break;

    case 'c':
      combined_loop_g_ = 1;
      simple_loop_g_ = 0;
      openmp_loop_g_ = 0;
      cache_loop_g_ = 0;
      combined_simple_mat_loop_g_ = 0;
      break;

    case 'm':
      combined_loop_g_ = 0;
      simple_loop_g_ = 0;
      openmp_loop_g_ = 0;
      cache_loop_g_ = 0;
      combined_simple_mat_loop_g_ = 1;
      break;

    case 'x':
      x_size_g_ = atoi((char *)optarg);
      break;

    case 'y':
      y_size_g_ = atoi((char *)optarg);
      break;

    case 'z':
      z_size_g_ = atoi((char *)optarg);
      break;

    case 'n':
      n_g_ = atoi((char *)optarg);
      break;
    }

  return 0;
}


/****************************************************************
 * MAIN
 ****************************************************************/
int main(int argc, char **argv)
{
  field_t eps[2], mu[2];
  unsigned int num_time_steps, i, j, k; 
  FILE *fields;
  unsigned int x, y, z, numthreads;
  int idx;
  time_t start, now;
  clock_t cpu_start, cpu_now;

  simple_loop_g_ = 1;
  openmp_loop_g_ = 0;
  cache_padding_g_ = 0;
  x_size_g_ = y_size_g_ = 128;
  z_size_g_ = 192;
  n_g_ = 100;
  quiet_g_ = 0;
  combined_loop_g_ = 0;

  decode_switches(argc, argv);

  x = 20;
  y = 20;
  z = 15;

  deltax_ = deltay_ = deltaz_ = 18.75e-9;
  deltat_ = 36e-18;

  dimx_ = x_size_g_; dimy_ = y_size_g_; dimz_ = z_size_g_;
  num_time_steps = n_g_;
  z_padding_ = 0;
  y_padding_ = 0;

  eps[0] = 1;
  eps[1] = 2.2;
  mu[0] = 1;
  mu[1] = 1;

  printf("High speed C FDTD test program\n");
  printf("Grid is %ix%ix%i, running for %i time steps.\n\n",
         dimx_, dimy_, dimz_, num_time_steps);
  
  load_materials(2, eps, mu);
  alloc_grid();

  for (i = 0; i < dimx_; i++)
    for (j = 0; j < dimy_; j++)
    {
      idx = pi(i,j,0);
      for (k = 0; k < dimz_; k++)
	material_[idx++] = 1;
    }

  /*fields = fopen("fields.txt", "w");*/

  if (openmp_loop_g_)
    printf("Using the OpenMP update loops.\n");

  if (simple_loop_g_)
    printf("Using the simple update loops.\n");
  
  if (cache_loop_g_)
    printf("Using the loops that hopefully exploit cache better.\n");

  if (combined_loop_g_)
    printf("Using the combined updates; 3 field components at once.\n");

  if (combined_simple_mat_loop_g_)
    printf("Using the combined updates; with homogeneuous material properties.\n");

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
    if (!quiet_g_)
      printf("High speed C, time step %i, source: %g\n", 
             i, gaussm(i, 500e12, 1, 300e12));

#ifdef USE_OPENMP
    if (openmp_loop_g_)
      omp_e_update();
#endif

    if (simple_loop_g_)
      e_update(); 

    if (cache_loop_g_)
      cache_e_update();

    if (combined_loop_g_)
      combined_e_update();

    if(combined_simple_mat_loop_g_)
      combined_simple_mat_e_update();

    ey_[pi(20, 20, 20)] = ey_[pi(20, 20, 20)] + gaussm(i, 500e12, 1, 300e12);

#ifdef USE_OPENMP
    if (openmp_loop_g_)
      omp_h_update();
#endif

    /* restricted_h_update(); */
    if (simple_loop_g_)
      h_update(); 

    if (cache_loop_g_)
      cache_h_update();

    if (combined_loop_g_)
      combined_h_update();

    if(combined_simple_mat_loop_g_)
      combined_simple_mat_h_update();

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
  unsigned int zbytes;
  uintptr_t e;
  int idx1, idx2, idx3;
  field_t *ez, *hy1, *hy2, *hx1, *hx2;

  sz = dimx_ * dimy_ * dimz_ * sizeof(field_t) + 16;

  temp = (field_t *)malloc(sizeof(field_t) * dimz_ + 16);

  /* Now with padding for cache considerations */ 
  if (cache_padding_g_)
  {
    sz = sz + 301; 
  }

  if (sz > 0) 
  {
    /* The allocated memory should be aligned on 128 byte boundaries,
       and each should be arranged so that it starts a different cache
       line than the rest. */ 

    /* Fortunatly for us, on OS X malloc returns memory block aligned
       on 16 byte boundaries required by AltiVec */ 
    ex_orig_ = (field_t *)malloc(sz);
    ey_orig_ = (field_t *)malloc(sz);
    ez_orig_ = (field_t *)malloc(sz);
    
    hx_orig_ = (field_t *)malloc(sz);
    hy_orig_ = (field_t *)malloc(sz);
    hz_orig_ = (field_t *)malloc(sz);
    
    printf("Allocated offset from 16 byte boundaries for origional ex, ey, ez, etc:\n");
    printf("0x%x, 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n",
           (uintptr_t)ex_orig_ & 0xf, 
           (uintptr_t)ey_orig_ & 0xf, 
           (uintptr_t)ez_orig_ & 0xf, 
           (uintptr_t)hx_orig_ & 0xf, 
           (uintptr_t)hy_orig_ & 0xf, 
           (uintptr_t)hz_orig_ & 0xf);

    zbytes = dimz_ * sizeof(mat_coef_t) + 16;
    Ca_temp_orig_ = (mat_coef_t *)malloc(zbytes);
    Cbx_temp_orig_ = (mat_coef_t *)malloc(zbytes);
    Cby_temp_orig_ = (mat_coef_t *)malloc(zbytes);
    Cbz_temp_orig_ = (mat_coef_t *)malloc(zbytes);

    Da_temp_orig_ = (mat_coef_t *)malloc(zbytes);
    Dbx_temp_orig_ = (mat_coef_t *)malloc(zbytes);
    Dby_temp_orig_ = (mat_coef_t *)malloc(zbytes);
    Dbz_temp_orig_ = (mat_coef_t *)malloc(zbytes);

    Ca_temp_ = Ca_temp_orig_ + ((uintptr_t)Ca_temp_orig_ & 0x0f) 
      / sizeof(mat_coef_t);
    Cbx_temp_ = Cbx_temp_orig_ + ((uintptr_t)Cbx_temp_orig_ & 0x0f) 
      / sizeof(mat_coef_t);
    Cby_temp_ = Cby_temp_orig_ + ((uintptr_t)Cby_temp_orig_ & 0x0f) 
      / sizeof(mat_coef_t);
    Cbz_temp_ = Cbz_temp_orig_ + ((uintptr_t)Cbz_temp_orig_ & 0x0f) 
      / sizeof(mat_coef_t);

    Da_temp_ = Da_temp_orig_ + ((uintptr_t)Da_temp_orig_ & 0x0f) 
      / sizeof(mat_coef_t);
    Dbx_temp_ = Dbx_temp_orig_ + ((uintptr_t)Dbx_temp_orig_ & 0x0f) 
      / sizeof(mat_coef_t);
    Dby_temp_ = Dby_temp_orig_ + ((uintptr_t)Dby_temp_orig_ & 0x0f) 
      / sizeof(mat_coef_t);
    Dbz_temp_ = Dbz_temp_orig_ + ((uintptr_t)Dbz_temp_orig_ & 0x0f) 
      / sizeof(mat_coef_t);

    printf("Allocated offset from 16 byte boundaries for origional Ca, Da, etc:\n");
    printf("0x%x, 0x%x, 0x%x, 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n",
           (uintptr_t)Ca_temp_ & 0xf, 
           (uintptr_t)Cbx_temp_ & 0xf, 
           (uintptr_t)Cby_temp_ & 0xf, 
           (uintptr_t)Cbz_temp_ & 0xf, 
           (uintptr_t)Da_temp_ & 0xf, 
           (uintptr_t)Dbx_temp_ & 0xf,
           (uintptr_t)Dby_temp_ & 0xf,
           (uintptr_t)Dbz_temp_ & 0xf);

    material_ = (unsigned int *)malloc(sz);

    if (ex_orig_ && ey_orig_ && ez_orig_ && hx_orig_ && hy_orig_ 
        && hz_orig_ && material_
        && Ca_temp_ && Cbx_temp_ && Cby_temp_ && Cbz_temp_
        && Da_temp_ && Dbx_temp_ && Dby_temp_ && Dbz_temp_)
    {
      memset(ex_orig_, 0, sz);
      memset(ey_orig_, 0, sz);
      memset(ez_orig_, 0, sz);
      
      memset(hx_orig_, 0, sz);
      memset(hy_orig_, 0, sz);
      memset(hz_orig_, 0, sz);

      memset(material_, 0, sz);

      memset(Ca_temp_orig_, 0, zbytes);
      memset(Cbx_temp_orig_, 0, zbytes);
      memset(Cby_temp_orig_, 0, zbytes);
      memset(Cbz_temp_orig_, 0, zbytes);

      memset(Da_temp_orig_, 0, zbytes);
      memset(Dbx_temp_orig_, 0, zbytes);
      memset(Dby_temp_orig_, 0, zbytes);
      memset(Dbz_temp_orig_, 0, zbytes);

    } else {
      printf("Insufficient memory; try reducing the grid size. \n");
      exit(1);
    }
    
    ex_ = ex_orig_ + ((uintptr_t)ex_orig_ & 0x0f) / sizeof(field_t);
    ey_ = ey_orig_ + ((uintptr_t)ey_orig_ & 0x0f) / sizeof(field_t);
    ez_ = ez_orig_ + ((uintptr_t)ez_orig_ & 0x0f) / sizeof(field_t);

    hx_ = hx_orig_ + ((uintptr_t)hx_orig_ & 0x0f) / sizeof(field_t);
    hy_ = hy_orig_ + ((uintptr_t)hy_orig_ & 0x0f) / sizeof(field_t);
    hz_ = hz_orig_ + ((uintptr_t)hz_orig_ & 0x0f) / sizeof(field_t);

    printf("Allocated offset from 16 byte boundaries for ex, ey, ez, etc:\n");
    printf("0x%x, 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n",
           (uintptr_t)ex_ & 0xf, 
           (uintptr_t)ey_ & 0xf, 
           (uintptr_t)ez_ & 0xf, 
           (uintptr_t)hx_ & 0xf, 
           (uintptr_t)hy_ & 0xf, 
           (uintptr_t)hz_ & 0xf);

    printf("No padding for cache:\n");
    e = (uintptr_t)ex_;
    printf("ex_orig ptr and low 22 bits of ex_orig pointer: %p, %p\n", 
           (void *)e, (void *)(e & 0x3FFFFF));
    e = (uintptr_t)ey_;
    printf("ey_orig ptr and low 22 bits of ey_orig pointer: %p, %p\n", 
           (void *)e, (void *)(e & 0x3FFFFF));
    e = (uintptr_t)ez_;
    printf("ez_orig ptr and low 22 bits of ez_orig pointer: %p, %p\n", 
           (void *)e, (void *)(e & 0x3FFFFF));

    e = (uintptr_t)hx_;
    printf("hx_orig ptr and low 22 bits of hx_orig pointer: %p, %p\n", 
           (void *)e, (void *)(e & 0x3FFFFF));
    e = (uintptr_t)hy_;
    printf("hy_orig ptr and low 22 bits of hy_orig pointer: %p, %p\n", 
           (void *)e, (void *)(e & 0x3FFFFF));
    e = (uintptr_t)hz_;
    printf("hz_orig ptr and low 22 bits of hz_orig pointer: %p, %p\n", 
           (void *)e, (void *)(e & 0x3FFFFF));
    
    idx1 = pi(16, 16, 0);
    idx2 = pi(15, 16, 0);
    idx3 = pi(16, 15, 0);
    
    ez = &ez_[idx1];
    hy1 = &hy_[idx1];
    hy2 = &hy_[idx2];
    hx1 = &hx_[idx3];
    hx2 = &hx_[idx1];

    e = (uintptr_t)ez;
    printf("ez ptr and low 22 bits of ez pointer: %p, %p\n", 
           (void *)e, (void *)(e & 0x3FFFFF));
    e = (uintptr_t)hy1;
    printf("hy1 ptr and low 22 bits of hy1 pointer: %p, %p\n", 
           (void *)e, (void *)(e & 0x3FFFFF));
    e = (uintptr_t)hy2;
    printf("hy2 ptr and low 22 bits of hy2 pointer: %p, %p\n", 
           (void *)e, (void *)(e & 0x3FFFFF));
    e = (uintptr_t)hx1;
    printf("hx1 ptr and low 22 bits of hx1 pointer: %p, %p\n", 
           (void *)e, (void *)(e & 0x3FFFFF));
    e = (uintptr_t)hx2;
    printf("hx2 ptr and low 22 bits of hx2 pointer: %p, %p\n", 
           (void *)e, (void *)(e & 0x3FFFFF));
  }
}

void free_grid()
{
  if (ex_ || ey_ || ez_ || hx_ || hy_ || hz_) {
    free(ex_orig_);
    free(ey_orig_);
    free(ez_orig_);

    free(hx_orig_);
    free(hy_orig_);
    free(hz_orig_);

    free(material_);

    free(Ca_temp_orig_);
    free(Cbx_temp_orig_);
    free(Cby_temp_orig_);
    free(Cbz_temp_orig_);

    free(Da_temp_orig_);
    free(Dbx_temp_orig_);
    free(Dby_temp_orig_);
    free(Dbz_temp_orig_);
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

  Ca_ = (mat_coef_t *)malloc(sz);
  Da_ = (mat_coef_t *)malloc(sz);

  Cbx_ = (mat_coef_t *)malloc(sz);
  Dbx_ = (mat_coef_t *)malloc(sz);

  /* Save some memory if possible. */
  if (deltay_ == deltax_)
  {
    Cby_ = Cbx_;
    Dby_ = Dbx_;
  } else {
    Cby_ = (mat_coef_t *)malloc(sz);
    Dby_ = (mat_coef_t *)malloc(sz);
  }

  if (deltaz_ == deltax_)
  {
    Cbz_ = Cbx_;
    Dbz_ = Dbx_;
  } else if (deltaz_ == deltay_) {
    Cbz_ = Cby_;
    Dbz_ = Dby_;
  } else {
    Cbz_ = (mat_coef_t *)malloc(sz);
    Dbz_ = (mat_coef_t *)malloc(sz);
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




