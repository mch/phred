/* Simple pointer arithmetic update equations using restricted pointer
   aliasing */

#include "common.h"

static void update_e_common();
/* static void update_ey(); */
/* static void update_ez(); */

static void update_h_common();
/* static void update_hy(); */
/* static void update_hz(); */

void restricted_e_update()
{
  /* update_ex(); */
/*   update_ey(); */
/*   update_ez(); */
  update_e_common();
}

void restricted_h_update()
{
/*   update_hx(); */
/*   update_hy(); */
/*   update_hz(); */
  update_h_common();
}

void update_e_common()
{
  unsigned int mid, idx, idxi, idxj, i, j, k;
  field_t * __restrict__ ex, * __restrict__ ey, * __restrict__ ez, 
    * __restrict__ hx, * __restrict__ hy, * __restrict__ hz, 
    * __restrict__ y_hz1, * __restrict__ x_hz2, * __restrict__ z_hy2, 
    * __restrict__ z_hx1;

  for (i = 1; i < dimx_ - 1; i++) {
    for (j = 1; j < dimy_ - 1; j++) {
      
      idx = pi(i, j, 1);
      idxj = pi(i, j-1, 1);
      idxi = pi(i-1, j, 1);

      ex = &(ex_[idx]);
      ey = &(ey_[idx]);
      ez = &(ez_[idx]);

      hz = &(hz_[idx]);
      x_hz2 = &(hz_[idxj]);
      y_hz1 = &(hz_[idxi]);
     
      hy = &(hy_[idx]);
      z_hy2 = &(hy_[pi(i-1, j, 0)]);

      hx = &(hx_[idx]);
      z_hx1 = &(hx_[pi(i, j-1, 0)]);

      // non-contiguous memory access' are consolodated here
      for (k = 0; k < dimz_ - 1; k++, idx++)
      {
        mid = material_[idx];

        Ca_temp_[k] = Ca_[mid];
        Cbx_temp_[k] = Cbx_[mid];
        Cby_temp_[k] = Cby_[mid];
        Cbz_temp_[k] = Cbz_[mid];
      }

      // This loop can now be vectorized
      for (k = 0; k < dimz_ - 1; k++) {
        ex[k] = Ca_temp_[k] * ex[k]
          + Cby_temp_[k] * (hz[k] - x_hz2[k])
          + Cbz_temp_[k] * (hy[k - 1] - hy[k]);

        ey[k] = Ca_temp_[k] * ey[k]
          + Cbz_temp_[k] * (hx[k] - hx[k-1])
          + Cbx_temp_[k] * (y_hz1[k] - hz[k]);          

        ez[k] = Ca_temp_[k] * ez[k]
          + Cbx_temp_[k] * (hy[k] - z_hy2[k])
          + Cby_temp_[k] * (z_hx1[k] - hx[k]);
      }
    }
  }

  // FIX THE EDGE CASES:
  
}


void update_h_common()
{
  unsigned int mid, idx, idxi, idxj, i, j, k;
  field_t * __restrict__ ex, * __restrict__ ey, * __restrict__ ez, 
    * __restrict__ hx, * __restrict__ hy, * __restrict__ hz, 
    * __restrict__ y_ez1, * __restrict__ x_ez2, * __restrict__ z_ey2, 
    * __restrict__ z_ex1;

  for (i = 1; i < dimx_ - 1; i++) {
    for (j = 1; j < dimy_ - 1; j++) {
      
      idx = pi(i, j, 0);
      idxi = pi(i+1, j, 0);
      idxj = pi(i, j+1, 0);
      
      hx = &(hx_[idx]);
      hy = &(hy_[idx]);
      hz = &(hz_[idx]);

      ex = &(ex_[idx]);
      z_ex1 = &(ex_[pi(i, j+1, 0)]);

      ey = &(ey_[idx]);
      z_ey2 = &(ey_[pi(i+1, j, 0)]);

      ez = &(ez_[idx]);
      x_ez2 = &(ez_[pi(i, j+1, 0)]);
      y_ez1 = &(ez_[pi(i+1, j, 0)]);

      // non-contiguous memory access' are consolodated here
      for (k = 0; k < dimz_ - 1; k++, idx++)
      {
        mid = material_[idx];

        Da_temp_[k] = Da_[mid];
        Dbx_temp_[k] = Dbx_[mid];
        Dby_temp_[k] = Dby_[mid];
        Dbz_temp_[k] = Dbz_[mid];
      }

      // This loop can now be vectorized
      for (k = 0; k < dimz_ - 1; k++) {
        hx[k] = Da_temp_[k] * hx[k]
          + Dby_temp_[k] * (ez[k] - x_ez2[k])
          + Dbz_temp_[k] * (ey[k+1] - ey[k]);

        hy[k] = Da_temp_[k] * hy[k]
          + Dbz_temp_[k] * (ex[k] - ex[k + 1])
          + Dbx_temp_[k] * (y_ez1[k] - ez[k]);           

        hz[k] = Da_temp_[k] * hz[k]
          + Dbx_temp_[k] * (ey[k] - z_ey2[k])
          + Dby_temp_[k] * (z_ex1[k] - ex[k]);
      }
    }
  }

  // FIX THE EDGE CASES:

}

/* void update_ex()  */
/* { */
/*   unsigned int mid, idx, idx2, i, j, k; */
/*   field_t * __restrict__ ex, * __restrict__ hz1, * __restrict__ hz2, * __restrict__ hy; */

/*   for (i = 0; i < dimx_; i++) { */
/*     for (j = 1; j < dimy_; j++) { */
      
/*       idx = pi(i, j, 1); */
/*       idx2 = pi(i, j-1, 1); */
/*       ex = &(ex_[idx]); */
/*       hz1 = &(hz_[idx]); */
/*       hz2 = &(hz_[idx2]); */
/*       hy = &(hy_[idx]); */

/*       for (k = 0; k < dimz_ - 1; k++, idx++) */
/*       { */
/*         mid = material_[idx]; */
/*         Coef0[k] = Ca_[mid]; */
/*         Coef1[k] = Cby_[mid]; */
/*         Coef2[k] = Cbz_[mid]; */
/*       } */

/*       for (k = 0; k < dimz_ - 1; k++) { */
/*         ex[k] = Coef0[k] * ex[k] */
/*           + Coef1[k] * (hz1[k] - hz2[k]) */
/*           + Coef2[k] * (hy[k - 1] - hy[k]); */

/*         //ex++; */
/*         //hz1++; */
/*         //hz2++; */
/*         //hy++; */
/*         //idx++; */
/*       } */
/*     } */
/*   } */
/* } */

/* void update_ey()  */
/* { */
/*   unsigned int mid, i, j, k, idx; */
/*   field_t * __restrict__ ey, * __restrict__ hx, * __restrict__ hz1, * __restrict__ hz2; */

/*   for (i = 1; i < dimx_; i++) { */
/*     for (j = 0; j < dimy_; j++) { */

/*       idx = pi(i, j, 1); */
/*       ey = &(ey_[idx]); */
/*       hx = &(hx_[idx]); */
/*       hz1 = &(hz_[pi(i-1, j, 1)]); */
/*       hz2 = &(hz_[idx]); */

/*       for (k = 0; k < dimz_ - 1; k++, idx++) */
/*       { */
/*         mid = material_[idx]; */
/*         Coef0[k] = Ca_[mid]; */
/*         Coef1[k] = Cbz_[mid]; */
/*         Coef2[k] = Cbx_[mid]; */
/*       } */

/*       for (k = 0; k < dimz_ - 1; k++) { */
/*         ey[k] = Coef0[k] * ey[k] */
/*           + Coef1[k] * (hx[k] - hx[k-1]) */
/*           + Coef2[k] * (hz1[k] - hz2[k]); */
/*       } */
/*     } */
/*   } */

/* } */

/* void update_ez()  */
/* { */
/*   unsigned int mid, i, j, k, idx; */
/*   field_t * __restrict__ ez, * __restrict__ hy1, * __restrict__ hy2,  */
/*     * __restrict__ hx1, * __restrict__ hx2; */
  
/*   for (i = 1; i < dimx_; i++) { */
/*     for (j = 1; j < dimy_; j++) { */

/*       idx = pi(i, j, 0); */
/*       ez = &(ez_[idx]); */
/*       hy1 = &(hy_[idx]); */
/*       hy2 = &(hy_[pi(i-1, j, 0)]); */
/*       hx1 = &(hx_[pi(i, j-1, 0)]); */
/*       hx2 = &(hx_[idx]); */

/*       for (k = 0; k < dimz_; k++, idx++) */
/*       { */
/*         mid = material_[idx]; */
/*         Coef0[k] = Ca_[mid]; */
/*         Coef1[k] = Cbx_[mid]; */
/*         Coef2[k] = Cby_[mid]; */
/*       } */

/*       for (k = 0; k < dimz_; k++) { */
/*         ez[k] = Coef0[k] * ez[k] */
/*           + Coef1[k] * (hy1[k] - hy2[k]) */
/*           + Coef2[k] * (hx1[k] - hx2[k]); */
/*       } */
/*     } */
/*   } */
/* } */

/* void update_hx() */
/* { */
/*   unsigned int mid, i, j, k, idx; */
/*   field_t * __restrict__ hx, * __restrict__ ez1, * __restrict__ ez2, * __restrict__ ey; */

/*   for (i = 0; i < dimx_; i++) { */
/*     for (j = 0; j < dimy_ - 1; j++) { */

/*       idx = pi(i, j, 0); */
/*       hx = &(hx_[idx]); */
/*       ez1 = &(ez_[idx]); */
/*       ez2 = &(ez_[pi(i, j+1, 0)]); */
/*       ey = &(ey_[idx]); */

/*       for (k = 0; k < dimz_ - 1; k++, idx++) */
/*       { */
/*         mid = material_[idx]; */
/*         Coef0[k] = Da_[mid]; */
/*         Coef1[k] = Dby_[mid]; */
/*         Coef2[k] = Dbz_[mid]; */
/*       } */

/*       for (k = 0; k < dimz_ - 1; k++) { */
/*         hx[k] = Coef0[k] * hx[k] */
/*           + Coef1[k] * (ez1[k] - ez2[k]) */
/*           + Coef2[k] * (ey[k+1] - ey[k]); */
/*       } */
/*     } */
/*   } */
/* } */

/* void update_hy() */
/* { */
/*   unsigned int mid, i, j, k, idx; */
/*   field_t * __restrict__ hy, * __restrict__ ex, * __restrict__ ez1, * __restrict__ ez2; */

/*   for (i = 0; i < dimx_ - 1; i++) { */
/*     for (j = 0; j < dimy_; j++) { */

/*       idx = pi(i, j, 0); */
/*       hy = &(hy_[idx]); */
/*       ex = &(ex_[idx]); */
/*       ez1 = &(ez_[pi(i+1, j, 0)]); */
/*       ez2 = &(ez_[idx]); */

/*       for (k = 0; k < dimz_ - 1; k++, idx++) */
/*       { */
/*         mid = material_[idx]; */
/*         Coef0[k] = Da_[mid]; */
/*         Coef1[k] = Dbz_[mid]; */
/*         Coef2[k] = Dbx_[mid]; */
/*       } */

/*       for (k = 0; k < dimz_ - 1; k++) { */
/*         hy[k] = Coef0[k] * hy[k] */
/*           + Coef1[k] * (ex[k] - ex[k + 1]) */
/*           + Coef2[k] * (ez1[k] - ez2[k]);         */
/*       } */
/*     } */
/*   } */
/* } */

/* void update_hz() */
/* { */
/*   unsigned int mid, i, j, k, idx; */
/*   field_t * __restrict__ hz1, * __restrict__ ey1, * __restrict__ ey2,  */
/*     * __restrict__ ex1, * __restrict__ ex2; */

/*   for (i = 0; i < dimx_ - 1; i++) { */
/*     for (j = 0; j < dimy_ - 1; j++) { */

/*       idx = pi(i, j, 0); */
/*       hz1 = &(hz_[idx]); */
/*       ey1 = &(ey_[idx]); */
/*       ey2 = &(ey_[pi(i+1, j, 0)]); */
/*       ex1 = &(ex_[pi(i, j+1, 0)]); */
/*       ex2 = &(ex_[idx]); */

/*       for (k = 0; k < dimz_; k++, idx++) */
/*       { */
/*         mid = material_[idx]; */
/*         Coef0[k] = Da_[mid]; */
/*         Coef1[k] = Dbx_[mid]; */
/*         Coef2[k] = Dby_[mid]; */
/*       } */

/*       for (k = 0; k < dimz_; k++) { */
/*         hz1[k] = Coef0[k] * hz1[k] */
/*           + Coef1[k] * (ey1[k] - ey2[k]) */
/*           + Coef2[k] * (ex1[k] - ex2[k]); */
/*       } */
/*     } */
/*   } */
/* } */

