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

unsigned int pi(unsigned int x, unsigned int y, 
                unsigned int z);

