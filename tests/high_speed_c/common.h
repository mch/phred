/* Common stuff needed by all update equation implementations */

#include "../../src/config.h"
#include "../../src/Types.hh"
#include "../../src/Constants.hh"

/****************************************************************
 * Global Variables!
 *
 * Not only are globals bad style, they can impede the performace of
 * high performance algorithms. A compiler can't keep a global
 * variable in a register the way it can with a local variable,
 * because of the possibility that the value may have been changed by
 * an interrupt handler or by another thread. 
 *
 * We're ok here though, because these pointers are memory bound and
 * always have to be loaded in anyway. 
 ****************************************************************/

/* Number of materials we know about (0 is PEC)*/
extern unsigned int num_materials_;

/* E Field Material Coefficients */
extern mat_coef_t *Ca_;
extern mat_coef_t *Cbx_;
extern mat_coef_t *Cby_;
extern mat_coef_t *Cbz_;

/* H Field Coefficients */
extern mat_coef_t *Da_;
extern mat_coef_t *Dbx_;
extern mat_coef_t *Dby_;
extern mat_coef_t *Dbz_;

/* Temporary coefficient holders. Pulls all the non-contiguous memory
   access out of the main update loop... hopefully making it faster. */
extern mat_coef_t * restrict Ca_temp_;
extern mat_coef_t * restrict Cbx_temp_;
extern mat_coef_t * restrict Cby_temp_;
extern mat_coef_t * restrict Cbz_temp_;

extern mat_coef_t * restrict Da_temp_;
extern mat_coef_t * restrict Dbx_temp_;
extern mat_coef_t * restrict Dby_temp_;
extern mat_coef_t * restrict Dbz_temp_;

/* The fields that get operated on. These are padded a bit so that
   better use of cache can be made. */ 
extern field_t *ex_;
extern field_t *ey_;
extern field_t *ez_;
extern field_t *hx_;
extern field_t *hy_;
extern field_t *hz_;

/* These are the origionally allocated pointers for field data. Used
   to free() the memory. */ 
extern field_t *ex_orig_;
extern field_t *ey_orig_;
extern field_t *ez_orig_;
extern field_t *hx_orig_;
extern field_t *hy_orig_;
extern field_t *hz_orig_;

/* Temporary data if needed */ 
extern field_t *temp;

/* The material for each point in the grid. This is an index into
 * the material arrays, Ca, Cbx, etc. */
extern unsigned int *material_;

/* Time and space steppings; the distance between each point in the
 * grid. */ 
extern delta_t deltax_;
extern delta_t deltay_;
extern delta_t deltaz_;
extern delta_t deltat_;

/* Size of the grid along each dimension */
extern unsigned int dimx_;
extern unsigned int dimy_;
extern unsigned int dimz_;


/****************************************************************
 * Option flags from the command line
 ****************************************************************/
extern char cache_padding_g_;


/****************************************************************
 * Function declarations
 ****************************************************************/

#ifdef USE_ALTIVEC
void av_e_update();
void av_h_update();
#endif

#ifdef USE_OPENMP
void omp_e_update();
void omp_h_update();

void omp_e_update2();
void omp_h_update2();
#endif

void e_update();
void h_update();

void restricted_e_update();
void restricted_h_update();

void cache_e_update();
void cache_h_update();

static inline unsigned int pi(unsigned int x, unsigned int y, 
                       unsigned int z)
{
  return z + (y + x*dimy_) * dimz_;
}

