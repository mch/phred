/* Simple pointer arithmetic update equations */

#include "common.h"

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
  field_t *ex, *hz1, *hz2, *hy, *temp2;

  for (i = 0; i < dimx_; i++) {
    for (j = 1; j < dimy_; j++) {
      
      idx = pi(i, j, 1);
      idx2 = pi(i, j-1, 1);
      ex = &(ex_[idx]);
      hz1 = &(hz_[idx]);
      hz2 = &(hz_[idx2]);
      hy = &(hy_[idx]);
      temp2 = temp;

      for (k = 0; k < dimz_ - 1; k++, idx++) {
        mid = material_[idx];
        
        *temp2 = /*Ca_[mid] */ ex[k]; // Hot, 12.83
        *temp2 += /*Cby_[mid] */ (hz1[k] - hz2[k]);
        *temp2 += /*Cbz_[mid] */ (hy[k - 1] - hy[k]);

        //ex++;
        //hz1++;
        //hz2++;
        //hy++;
        //idx++;
        temp2++;
      }
      memcpy(&ex_[pi(i, j, 0)], temp, sizeof(field_t) * dimz_);
    }
  }
}

void update_ey() 
{
  unsigned int mid, i, j, k, idx, *mat;
  field_t *ey, *hx, *hz1, *hz2, *temp2;

  for (i = 1; i < dimx_; i++) {
    for (j = 0; j < dimy_; j++) {

      idx = pi(i, j, 1);
      ey = &(ey_[idx]);
      hx = &(hx_[idx]);
      hz1 = &(hz_[pi(i-1, j, 1)]);
      hz2 = &(hz_[idx]);
      mat = &material_[idx];

      temp2 = temp;
      for (k = 1; k < dimz_; k++) {
        mid = *mat++;

        *temp2 = /*Ca_[mid] * */ *ey;
        *temp2 += /*Cbz_[mid] * */(*hx - *(hx-1));
        *temp2 += /*Cbx_[mid] * */(*hz1 - *hz2);

        ey++;
        hx++;
        hz1++;
        hz2++;
        idx++;
        temp2++;
      }
      memcpy(&ey_[pi(i, j, 0)], temp, sizeof(field_t) * dimz_);
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
        mid = material_[idx];

        temp[k] = Ca_[mid] * *ez
          + Cbx_[mid] * (*hy1 - *hy2)
          + Cby_[mid] * (*hx1 - *hx2);

        ez++;
        hy1++; hy2++; hx1++; hx2++;
        idx++;
      }
      memcpy(&ez_[pi(i, j, 0)], temp, sizeof(field_t) * dimz_);
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

        temp[k] = Da_[mid] * *hx
          + Dby_[mid] * (*ez1 - *ez2)
          + Dbz_[mid] * (*(ey+1) - *ey);

        hx++; idx++;
        ez1++; ez2++; ey++;
      }
      memcpy(&hx_[pi(i, j, 0)], temp, sizeof(field_t) * dimz_);
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

        temp[k] = Da_[mid] * *hy // Hot, 19.39
          + Dbz_[mid] * (*ex - *(ex + 1))
          + Dbx_[mid] * (*ez1 - *ez2);        
        
        hy++; idx++;
        ex++; ez1++; ez2++;
      }
      memcpy(&hy_[pi(i, j, 0)], temp, sizeof(field_t) * dimz_);
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

        temp[k] = Da_[mid] * *hz1
          + Dbx_[mid] * (*ey1 - *ey2)
          + Dby_[mid] * (*ex1 - *ex2);

        hz1++; idx++;
        ey1++; ey2++;
        ex1++; ex2++;
      }
      memcpy(&hz_[pi(i, j, 0)], temp, sizeof(field_t) * dimz_);
    }
  }
}

