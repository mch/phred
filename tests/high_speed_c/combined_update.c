/* Updated all three field components at the same time. Currently
   cheats a bit, since it doesn't update the edges like it would
   normally have to. */ 

#include "common.h"

#include <string.h>
#include <inttypes.h>

static void update_ex();
static void update_ey();
static void update_ez();

static void update_hx();
static void update_hy();
static void update_hz();

void combined_e_update()
{
  unsigned int mid;
  int i, j, k, idx, idx_x, idx_y, idx_z, low_bound, up_bound;
  field_t *ex, *ey, *ez;
  field_t *hx_y, *hx_z;
  field_t *hy_x, *hy_z;
  field_t *hz, *hz_x, *hz_y;

  float *ca;
  float *cbx;
  float *cby;
  float *cbz;

  // Calculate the upper and lower bounds of the iteration over the z
  // array such that the starting address is aligned on a 16 byte
  // boundary and the span is a multiple of 16. 
  low_bound = 4;
  up_bound = dimz_ - 1;
  up_bound = up_bound - ((up_bound - low_bound) % 4);

  for (i = 1; i < dimx_; i++) {
    for (j = 1; j < dimy_; j++) {
      
      idx = pi(i, j, low_bound);
      idx_x = pi(i-1, j, low_bound);
      idx_y = pi(i, j-1, low_bound);
      idx_z = pi(i, j, low_bound - 1);

      ex = &ex_[idx];
      ey = &ey_[idx];
      ez = &ez_[idx];
      
      __assume_aligned(ex, 16);
      __assume_aligned(ey, 16);
      __assume_aligned(ez, 16);

      hz = &hz_[idx];

      hy_x = &hy_[idx_x];
      hz_x = &hz_[idx_x];

      hx_y = &hx_[idx_y];
      hz_y = &hz_[idx_y];

      __assume_aligned(hz, 16);
      __assume_aligned(hy_x, 16);
      __assume_aligned(hz_x, 16);
      __assume_aligned(hx_y, 16);
      __assume_aligned(hz_y, 16);

      hx_z = &hx_[idx_z];
      hy_z = &hy_[idx_z];

      up_bound = dimz_ - 1;

      for (k = 1; k < up_bound; k++) {
        mid = material_[idx++];
        Ca_temp_[k] = Ca_[mid];
        Cbx_temp_[k] = Cbx_[mid];
        Cby_temp_[k] = Cby_[mid];
        Cbz_temp_[k] = Cbz_[mid];
      }

      ca = &Ca_temp_[low_bound];
      cbx = &Cbx_temp_[low_bound];
      cby = &Cby_temp_[low_bound];
      cbz = &Cbz_temp_[low_bound];

      #pragma ivdep
      for (k = low_bound; k < up_bound; k++) {
        *ex = *ca * *ex
          + *cby * (*hz - *hz_y)
          + *cbz * (*hy_z - *(hy_z + 1));

        ca++; cby++; cbz++;
        ex++;

        hz++; hz_y++;
        hy_z++; 
      }

      ca = Ca_temp_;
      cbx = Cbx_temp_;
      cby = Cby_temp_;
      cbz = Cbz_temp_;

      #pragma ivdep
      for (k = low_bound; k < up_bound; k++) {
        *ey = *ca * *ey
          + *cbz * (*(hx_z + 1) - *hx_z)
          + *cbx * (*hz_x - *hz);

        ca++; cbx++; cbz++;
        ey++;

        hz++; hz_x++;
        hx_z++;
      }

      ca = Ca_temp_;
      cbx = Cbx_temp_;
      cby = Cby_temp_;
      cbz = Cbz_temp_;

      #pragma ivdep
      for (k = low_bound; k < up_bound; k++) {
        *ez = *ca * *ez
          + *cbx * (*(hy_z + 1) - *hy_x) 
          + *cby * (*hx_y - *(hx_z + 1));

        ca++; cbx++; cby++;
        ez++;

        hy_z++; hy_x++;
        hx_z++; hx_y++;
      }
    }
  }
}

void combined_h_update()
{
  unsigned int mid;
  int i, j, k, idx, idx_x, idx_y, idx_z, low_bound, up_bound;
  field_t *hx, *hy, *hz;
  field_t *ex, *ex_y;
  field_t *ey, *ey_x;
  field_t *ez, *ez_x, *ez_y;

  float *da;
  float *dbx;
  float *dby;
  float *dbz;

  // Calculate the upper and lower bounds of the iteration over the z
  // array such that the starting address is aligned on a 16 byte
  // boundary and the span is a multiple of 16. 
  low_bound = 4;
  up_bound = dimz_ - 1;

  for (i = 1; i < dimx_; i++) {
    for (j = 1; j < dimy_; j++) {
      
      idx = pi(i, j, low_bound);
      idx_x = pi(i+1, j, low_bound);
      idx_y = pi(i, j+1, low_bound);
      idx_z = pi(i, j, low_bound + 1);

      hx = &hx_[idx];
      hy = &hy_[idx];
      hz = &hz_[idx];

      ex = &ex_[idx];
      ey = &ey_[idx];
      ez = &ez_[idx];

      ey_x = &ey_[idx_x];
      ez_x = &ez_[idx_x];

      ex_y = &ex_[idx_y];
      ez_y = &ez_[idx_y];

      for (k = low_bound; k < up_bound; k++) {
        mid = material_[idx++];
        Da_temp_[k] = Da_[mid];
        Dbx_temp_[k] = Dbx_[mid];
        Dby_temp_[k] = Dby_[mid];
        Dbz_temp_[k] = Dbz_[mid];
      }

      da = Da_temp_;
      dbx = Dbx_temp_;
      dby = Dby_temp_;
      dbz = Dbz_temp_;

      #pragma ivdep
      for (k = low_bound; k < up_bound; k++) {
        *hx = *da * *hx
          + *dby * (*ez - *ez_y)
          + *dbz * (*(ey + 1) - *ey);

        da++; dbx++; dby++; dbz++;
        hx++;

        ez++; ez_y++; ey++;
      }

      da = Da_temp_;
      dbx = Dbx_temp_;
      dby = Dby_temp_;
      dbz = Dbz_temp_;

      #pragma ivdep
      for (k = low_bound; k < up_bound; k++) {
        *hy = *da * *hy
          + *dbz * (*ex - *(ex + 1))
          + *dbx * (*ez_x - *ez);

        da++; dbx++; dby++; dbz++;
        hy++;

        ex++; ez_x++; ez++;
      }

      da = Da_temp_;
      dbx = Dbx_temp_;
      dby = Dby_temp_;
      dbz = Dbz_temp_;

      #pragma ivdep
      for (k = low_bound; k < up_bound; k++) {
        *hz = *da * *hz
          + *dbx * (*ey - *ey_x)
          + *dby * (*ex_y - *ex);

        da++; dbx++; dby++; dbz++;
        hz++;

        ex++; ey++;
        ey_x++;
        ex_y++;
      }
    }
  }
}
