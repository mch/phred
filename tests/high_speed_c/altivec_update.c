/* Altivec update equations */

#include "common.h"

#ifdef USE_ALTIVEC
/** AltiVec Implementations, G4's and G5's only **/

#include "../../src/AltiVec.hh"

/* Individual component updates */
static void av_update_ex();
static void av_update_ey();
static void av_update_ez();
static void av_update_hx();
static void av_update_hy();
static void av_update_hz();

void av_e_update()
{
  av_update_ex();
  av_update_ey();
  av_update_ez();
}

void av_h_update()
{
  av_update_hx();
  av_update_hy();
  av_update_hz();  
}

/**
 * Computes the E field components using Altivec. 
 *
 * Ideas for optimization: 
 * - have one function that computes all E components at once, and
 *   one that computes all H components at once. The region will have
 *   to be reduced so that we don't hit the edges, but corrections
 *   can be made for the leftovers after. 
 */
#ifdef DONTUSE
void av_e_update()
{
  unsigned int mid, idx, idx2, i, j, k;

  vector float *ex, *hz1, *hz2, *hy1, *hy2, *hx1, *hx2;

  FloatVector_t Ca, Cby, Cbz;
  
  /*
   * Field components have to be aligned properly on boundaries for
   * fast loads and stores with AltiVec
   */

  /* Add one 1 to xmin, ymin, zmin, since only certian field
   * components are calculated at those points, and they can be done
   * as a correction after. */
  update_r.xmin++;
  update_r.ymin++;
  update_r.zmin++;

  /* It should be possible to eliminate these outer loops and have the
   * entire region updated by one loop, since the memory region is
   * contigious. There will be errors at the edges, but these can can
   * be corrected for after. If the result is faster...*/
  for (i = update_r.xmin; i < update_r.xmax; i++) {
    for (j = update_r.ymin; j < update_r.ymax; j++) {

      /* Start cache streams from memory. These computations are
         almost certianlly limited by memory, rather than by CPU. */ 
      int prefetch_const = 0x10010100; // 256 bytes
      /*vec_dst(ex, prefect_const, 0);*/ /* Start cache stream */

      idx = pi(i, j, update_r.zmin);
      idx2 = pi(i, j-1, update_r.zmin);
      ex = (vector float *)&(ex_[idx]);

      hx1 = (vector float *)&(hx_[pi(i, j-1, 0)]);
      hx2 = (vector float *)&(hx_[idx]);

      hy1 = (vector float *)&(hy_[idx]);
      hy2 = (vector float *)&(hy_[pi(i-1, j, 0)]);

      hz1 = (vector float *)&(hz_[idx]);
      hz2 = (vector float *)&(hz_[idx2]);

      /* Ensure we don't step over the edge of our region along z */
      unsigned int kmax = update_r.zmax - 4;
      for (k = update_r.zmin; k < kmax; k+=4) {
        /* Load material data into vectors (POTENTIAL MAJOR BOTTLENECK) */
        for (int l = 0; l < 4; l++)
        {
          mid = material_[idx + l];
          Ca.elements[l] = Ca_[mid];
          Cby.elements[l] = Cby_[mid];
          Cbz.elements[l] = Cbz_[mid];
        }

        /* Perform the actual computation */
        ex = vec_madd(Ca, ex, vec_madd(Cby, vec_sub(hz1, hz2), zero));
        ex = vec_madd(Cbz, vec_madd(hy1, hy2, zero), ex);

        ey = vec_madd(Ca, ey, vec_madd(Cby, vec_sub(hx1, hx2), zero));
        ey = vec_madd(Cbz, vec_madd(hz1, hz2, zero), ey);

        ez = vec_madd(Ca, ez, vec_madd(Cby, vec_sub(hy1, hy2), zero));
        ez = vec_madd(Cbz, vec_madd(hx1, hx2, zero), ez);

/*         *ex = Ca_[mid] * *ex */
/*           + Cby_[mid] * (*hz1 - *hz2) */
/*           + Cbz_[mid] * (*(hy - 1) - *hy); */

/*         *ey = Ca_[mid] * *ey */
/*           + Cbz_[mid] * (*hx - *(hx-1)) */
/*           + Cbx_[mid] * (*hz1 - *hz2); */

/*         *ez = Ca_[mid] * *ez */
/*           + Cbx_[mid] * (*hy1 - *hy2) */
/*           + Cby_[mid] * (*hx1 - *hx2); */

        ex++;
        hx1++;
        hx2++;
        hy1++;
        hy2++;
        hz1++;
        hz2++;
        idx += 4;
      }
      /*vec_dss(0);*/ /* Stop cache streaam */

      /* Do the 3 or fewer elements that can't be done in a vector */
/*       for (k = kmax; k < update_r.zmax; k++) { */
/*         mid = material_[idx];         */

/*         *ex = Ca_[mid] * *ex */
/*           + Cby_[mid] * (*hz1 - *hz2) */
/*           + Cbz_[mid] * (*(hy - 1) - *hy); */

/*         *ey = Ca_[mid] * *ey */
/*           + Cbz_[mid] * (*hx - *(hx-1)) */
/*           + Cbx_[mid] * (*hz1 - *hz2); */

/*         *ez = Ca_[mid] * *ez */
/*           + Cbx_[mid] * (*hy1 - *hy2) */
/*           + Cby_[mid] * (*hx1 - *hx2); */

/*         ex++; */
/*         hz1++; */
/*         hz2++; */
/*         hy++; */
/*         idx++; */
/*       } */
    }
  }
}
#endif

void av_update_ex() 
{
  unsigned int mid, idx, idx2, i, j, k;
  vector float *ex, *hz1, *hz2, *hy1, *hy2;
  FloatVector_t Ca, Cby, Cbz;
  unsigned int kmax;
  int prefetch_const = 0x10010100; // 256 bytes
  int l;

  for (i = 0; i < dimx_; i++) {
    for (j = 1; j < dimy_; j++) {
      
      /* HAVE TO MAKE SURE THE IDX's give pointers that are 16 byte
       * aligned. Corrections for data that is not on 16 byte
       * boundaries will have to be done separately.
       */
      idx = pi(i, j, 1);
      idx2 = pi(i, j-1, 1);
      ex = (vector float *)&(ex_[idx]);
      hz1 = (vector float *)&(hz_[idx]);
      hz2 = (vector float *)&(hz_[idx2]);
      hy = (vector float *)&(hy_[idx]);
      hy = (vector float *)&(hy_[idx]);

      /* Start cache streams from memory. These computations are
         almost certianlly limited by memory, rather than by CPU. */ 
      /*vec_dst(ex, prefect_const, 0);*/ /* Start cache stream */

      kmax = dimz_ - 4;
      for (k = 1; k < kmax; k+=4) {
        /* Load material data into vectors (POTENTIAL MAJOR BOTTLENECK) */
        for (l = 0; l < 4; l++)
        {
          mid = material_[idx + l];
          Ca.elements[l] = Ca_[mid];
          Cby.elements[l] = Cby_[mid];
          Cbz.elements[l] = Cbz_[mid];
        }

        /* AltiVec computation */
        ex = vec_madd(Ca, ex, vec_madd(Cby, vec_sub(hz1, hz2), ZERO));
        ex = vec_madd(Cbz, vec_madd(hy1, hy2, ZERO), ex);

        ex++;
        hz1++;
        hz2++;
        hy++;
        idx += 4;
      }

      /* Do the last 3 or fewer elements */
    }
  }
}

/* Straight out of Taflove. */
void av_update_ey() 
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
void av_update_ez() 
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

void av_update_hx()
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

void av_update_hy()
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

void av_update_hz()
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

#endif
