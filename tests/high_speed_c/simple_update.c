/* Simple pointer arithmetic update equations */

#include "common.h"

static void update_ex();
static void update_ey();
static void update_ez();

static void update_hx();
static void update_hy();
static void update_hz();

void e_update()
{
  update_ex();
  update_ey();
  update_ez();
}

void h_update()
{
  update_hx();
  update_hy();
  update_hz();
}

void update_ex() 
{
  unsigned int mid, idx, idx2, i, j, k;
  field_t *ex, *hz1, *hz2, *hy;

  for (i = 0; i < dimx_; i++) {
    for (j = 1; j < dimy_; j++) {
      
      idx = pi(i, j, 1);
      idx2 = pi(i, j-1, 1);
      ex = &(ex_[idx]);
      hz1 = &(hz_[idx]);
      hz2 = &(hz_[idx2]);
      hy = &(hy_[idx]);

      for (k = 1; k < dimz_; k++) {
        mid = material_[idx];
        
        *ex = Ca_[mid] * *ex
          + Cby_[mid] * (*hz1 - *hz2)
          + Cbz_[mid] * (*(hy - 1) - *hy);

        ex++;
        hz1++;
        hz2++;
        hy++;
        idx++;
      }
    }
  }
}

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

