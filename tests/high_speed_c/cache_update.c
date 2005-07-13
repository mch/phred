/* Simple pointer arithmetic update equations */

#include "common.h"

#include <string.h>

static void update_ex();
static void update_ey();
static void update_ez();

static void update_hx();
static void update_hy();
static void update_hz();

void cache_e_update()
{
  update_ex();
  update_ey();
  update_ez();
}

void cache_h_update()
{
  update_hx();
  update_hy();
  update_hz();
}

void update_ex() 
{
  unsigned int mid, idx, idx2;
  int i, j, k;
  field_t *ex, *hz1, *hz2, *hy, *hy2;

  for (i = 0; i < dimx_; i++) {
    for (j = 1; j < dimy_; j++) {
      
      idx = pi(i, j, 1);
      idx2 = pi(i, j-1, 1);
      ex = &(ex_[idx]);
      hz1 = &(hz_[idx]);
      hz2 = &(hz_[idx2]);
      hy = &(hy_[idx]);
      hy2 = hy - 1;

      int up_bound = dimz_ - 1;

      for (k = 0; k < up_bound; k++) {
        mid = material_[idx];
        Ca_temp_[k] = Ca_[mid];
        Cby_temp_[k] = Cby_[mid];
        Cbz_temp_[k] = Cbz_[mid];
      }

      float *ca = Ca_temp_;
      float *cby = Cby_temp_;
      float *cbz = Cbz_temp_;

      #pragma ivdep
      for (k = 0; k < up_bound; k++) {
        *ex = *ca * *ex // Hot, 12.83
          + *cby * (*hz1 - *hz2)
          + *cbz * (*hy2 - *hy);

        ca++;
        cby++;
        cbz++;
        ex++;
        hz1++;
        hz2++;
        hy++;
        hy2++;
        //temp2++;
      }
      //memcpy(&ex_[pi(i, j, 0)], temp, sizeof(field_t) * dimz_);
    }
  }
}

void update_ey() 
{
  unsigned int mid, i, j, k, idx;
  field_t *ey, *hx, *hx2, *hz1, *hz2;

  for (i = 1; i < dimx_; i++) {
    for (j = 0; j < dimy_; j++) {

      idx = pi(i, j, 1);
      ey = &(ey_[idx]);
      hx = &(hx_[idx]);
      hx2 = hx - 1;
      hz1 = &(hz_[pi(i-1, j, 1)]);
      hz2 = &(hz_[idx]);

      int up_bound = dimz_;

      for (k = 0; k < up_bound; k++) {
        mid = material_[idx];
        Ca_temp_[k] = Ca_[mid];
        Cbz_temp_[k] = Cbz_[mid];
        Cbx_temp_[k] = Cbx_[mid];
      }
      float *ca = Ca_temp_;
      float *cbz = Cbz_temp_;
      float *cbx = Cbx_temp_;

      #pragma ivdep
      for (k = 1; k < up_bound; k++) {
        *ey = *ca * *ey
          + *cbz * (*hx - *hx2)
          + *cbx * (*hz1 - *hz2);

        ey++; ca++; cbz++; cbz++;
        hx++;
        hz1++;
        hz2++;
        //temp2++;
      }
      //memcpy(&ey_[pi(i, j, 0)], temp, sizeof(field_t) * dimz_);
    }
  }

}

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

      int up_bound = dimz_;

      for (k = 0; k < up_bound; k++) {
        mid = material_[idx];
        Ca_temp_[k] = Ca_[mid];
        Cbx_temp_[k] = Cbx_[mid];
        Cby_temp_[k] = Cby_[mid];
      }
      float *ca = Ca_temp_;
      float *cby = Cby_temp_;
      float *cbx = Cbx_temp_;

      #pragma ivdep
      for (k = 0; k < up_bound; k++) {
        *ez = *ca * *ez
          + *cbx * (*hy1 - *hy2)
          + *cby * (*hx1 - *hx2);

        ca++; cbx++; cby++;
        ez++;
        hy1++; hy2++; hx1++; hx2++;
      }
      //memcpy(&ez_[pi(i, j, 0)], temp, sizeof(field_t) * dimz_);
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

      int up_bound = dimz_ - 1;

      for (k = 0; k < up_bound; k++) {
        mid = material_[idx];
        Da_temp_[k] = Da_[mid];
        Dbz_temp_[k] = Dbz_[mid];
        Dby_temp_[k] = Dby_[mid];
      }
      float *da = Da_temp_;
      float *dby = Dby_temp_;
      float *dbz = Dbz_temp_;

      #pragma ivdep
      for (k = 0; k < up_bound; k++) {
        *hx = *da * *hx
          + *dby * (*ez1 - *ez2)
          + *dbz * (*(ey+1) - *ey);

        da++; dby++; dbz++;
        hx++; idx++;
        ez1++; ez2++; ey++;
      }
      //memcpy(&hx_[pi(i, j, 0)], temp, sizeof(field_t) * dimz_);
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

      int up_bound = dimz_ - 1;

      for (k = 0; k < up_bound; k++) {
        mid = material_[idx];
        Da_temp_[k] = Da_[mid];
        Dbz_temp_[k] = Dbz_[mid];
        Dbx_temp_[k] = Dbx_[mid];
      }
      float *da = Da_temp_;
      float *dbx = Dbx_temp_;
      float *dbz = Dbz_temp_;

      #pragma ivdep
      for (k = 0; k < up_bound; k++) {
        *hy = *da * *hy // Hot, 19.39
          + *dbz * (*ex - *(ex + 1))
          + *dbx * (*ez1 - *ez2);        
        
        da++; dbz++; dbx++;
        hy++; idx++;
        ex++; ez1++; ez2++;
      }
      //memcpy(&hy_[pi(i, j, 0)], temp, sizeof(field_t) * dimz_);
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

      int up_bound = dimz_;

      for (k = 0; k < up_bound; k++) {
        mid = material_[idx];
        Da_temp_[k] = Da_[mid];
        Dbx_temp_[k] = Dbx_[mid];
        Dby_temp_[k] = Dby_[mid];
      }
      float *da = Da_temp_;
      float *dby = Dby_temp_;
      float *dbx = Dbx_temp_;

      #pragma ivdep
      for (k = 0; k < up_bound; k++) {
        *hz1 = *da * *hz1
          + *dbx * (*ey1 - *ey2)
          + *dby * (*ex1 - *ex2);

        da++; dbx++; dby++;
        hz1++; idx++;
        ey1++; ey2++;
        ex1++; ex2++;
      }
      //memcpy(&hz_[pi(i, j, 0)], temp, sizeof(field_t) * dimz_);
    }
  }
}

