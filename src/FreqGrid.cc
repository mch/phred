/* 
   phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004 Matt Hughes <mhughe@uvic.ca>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#include "config.h"

#include "FreqGrid.hh"
#include "Exceptions.hh"

#include <string.h> // for memset()
#include <math.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

FreqGrid::FreqGrid()
  : vcdt_(0), omegapsq_(0), dx_(0), sx_(0), sxm1_(0), sxm2_(0),
    dy_(0), sy_(0), sym1_(0), sym2_(0),
    dz_(0), sz_(0), szm1_(0), szm2_(0)
{}

FreqGrid::~FreqGrid()
{
}

void FreqGrid::alloc_grid()
{
  Grid::alloc_grid();

  // Ugly-be-gone: Now that the grid knows about the geometries, it
  // can check each and figure out which have to have aux data stored
  // in them and set them accordingly. 

  unsigned int sz = 0;

  sz = info_.dimx_ * info_.dimy_ * info_.dimz_;

  // UGLY!
  unsigned int num_plasma = 0;

  // TEMPORARY
  num_plasma = sz;

  if (num_plasma > 0)
  {
    dx_ = new field_t[num_plasma];
    sx_ = new field_t[num_plasma];
    sxm1_ = new field_t[num_plasma];
    sxm2_ = new field_t[num_plasma];
    
    dy_ = new field_t[num_plasma];
    sy_ = new field_t[num_plasma];
    sym1_ = new field_t[num_plasma];
    sym2_ = new field_t[num_plasma];
    
    dz_ = new field_t[num_plasma];
    sz_ = new field_t[num_plasma];
    szm1_ = new field_t[num_plasma];
    szm2_ = new field_t[num_plasma];
    
    if (!dx_ || !sx_ || !sxm1_ || !sxm2_
        || !dy_ || !sy_ || !sym1_ || !sym2_
        || !dz_ || !sz_ || !szm1_ || !szm2_)
    {
      free_grid();
      throw MemoryException();
    }

    memset(dx_, 0, sizeof(field_t) * num_plasma);
    memset(sx_, 0, sizeof(field_t) * num_plasma);
    memset(sxm1_, 0, sizeof(field_t) * num_plasma);
    memset(sxm2_, 0, sizeof(field_t) * num_plasma);

    memset(dy_, 0, sizeof(field_t) * num_plasma);
    memset(sy_, 0, sizeof(field_t) * num_plasma);
    memset(sym1_, 0, sizeof(field_t) * num_plasma);
    memset(sym2_, 0, sizeof(field_t) * num_plasma);

    memset(dz_, 0, sizeof(field_t) * num_plasma);
    memset(sz_, 0, sizeof(field_t) * num_plasma);
    memset(szm1_, 0, sizeof(field_t) * num_plasma);
    memset(szm2_, 0, sizeof(field_t) * num_plasma);
  }
}

void FreqGrid::free_grid()
{
  Grid::free_grid();

  if (dx_)
    delete[] dx_;
  
  if (sx_)
    delete[] sx_;
  
  if (sxm1_)
    delete[] sxm1_;

  if (sxm2_)
    delete[] sxm2_;

  if (dy_)
    delete[] dy_;
  
  if (sy_)
    delete[] sy_;
  
  if (sym1_)
    delete[] sym1_;

  if (sym2_)
    delete[] sym2_;

  if (dz_)
    delete[] dz_;
  
  if (sz_)
    delete[] sz_;
  
  if (szm1_)
    delete[] szm1_;

  if (szm2_)
    delete[] szm2_;

}
  
void FreqGrid::load_materials(MaterialLib &matlib)
{
  Grid::load_materials(matlib);

  int num_mat = matlib.num_materials() + 1;  
  vcdt_ = new mat_coef_t[num_mat];
  omegapsq_ = new mat_coef_t[num_mat];

  if (!vcdt_ || !omegapsq_) {
    free_material();
    throw MemoryException();
  }

  memset(vcdt_, 0, sizeof(mat_coef_t) * num_mat);
  memset(omegapsq_, 0, sizeof(mat_coef_t) * num_mat);

  int index = 0;

  vector<Material>::const_iterator iter = matlib.get_material_iter_begin();
  vector<Material>::const_iterator iter_e = matlib.get_material_iter_end();

  vcdt_[0] = 0;
  omegapsq_[0] = 0;
  ++index;

  while(iter != iter_e)
  {
    if ((*iter).get_collision_freq() > 0)
    {
      vcdt_[index] = exp(-1.0 * (*iter).get_collision_freq() * get_deltat());
      omegapsq_[index] = pow((*iter).get_plasma_freq(), 
                             static_cast<float>(2.0))
        * (get_deltat() / (*iter).get_collision_freq());
    } else {
      vcdt_[index] = 0.0;
      omegapsq_[index] = 0.0;
    }
    ++index;
    ++iter;
  }
}

void FreqGrid::free_material()
{
  Grid::free_material();

  if (vcdt_)
  {
    delete[] vcdt_;
    vcdt_ = 0;
  }

  if (omegapsq_)
  {
    delete[] omegapsq_;
    omegapsq_ = 0;
  }
}

void FreqGrid::update_ex(region_t update_r)
{
  unsigned int mid, idx, idx2, plasma_idx = 0;
  int i, j, k;
  field_t *ex, *hz1, *hz2, *hy;

  // Inner part
#ifdef USE_OPENMP
#pragma omp parallel private(mid, i, j, k, idx, idx2, ex, hz1, hz2, hy)
#endif
  {
#ifdef USE_OPENMP
#pragma omp for
#endif
    for (i = update_r.xmin; i < update_r.xmax; i++) {
      for (j = update_r.ymin; j < update_r.ymax; j++) {
      
        idx = pi(i, j, update_r.zmin);
        idx2 = pi(i, j-1, update_r.zmin);

        ex = &(ex_[idx]);
        hz1 = &(hz_[idx]);
        hz2 = &(hz_[idx2]);
        hy = &(hy_[idx]);
        
        for (k = update_r.zmin; k < update_r.zmax; k++) {
          mid = material_[idx];

          // Is this a plasma?
          if (vcdt_[mid] == 0.0) 
          {        
            *ex = Ca_[mid] * *ex
              + Cby_[mid] * (*hz1 - *hz2)
              + Cbz_[mid] * (*(hy - 1) - *hy);
          }
          else
          {
            *ex = Ca_[mid] * dx_[idx]
              + Cby_[mid] * (*hz1 - *hz2)
              + Cbz_[mid] * (*(hy - 1) - *hy);
            
            dx_[idx] = *ex;
            
            *ex = dx_[idx] - sx_[idx];
            
            sx_[idx] = (1 + vcdt_[mid]) * sxm1_[idx]
              - (vcdt_[mid] * sxm2_[idx])
              + (omegapsq_[mid] * (1 - vcdt_[mid])) * *ex;
            
            sxm2_[idx] = sxm1_[idx];
            sxm1_[idx] = sx_[idx];

            //++plasma_idx;
          }
          
          ex++;
          hz1++;
          hz2++;
          hy++;
          idx++;
        }
      }
    }
  }
}

void FreqGrid::update_ey(region_t update_r)
{
  unsigned int mid, idx, plasma_idx = 0;
  int i, j, k;
  field_t *ey, *hx, *hz1, *hz2;

  // Inner part
#ifdef USE_OPENMP
#pragma omp parallel private(mid, i, j, k, idx, ey, hx, hz1, hz2)
#endif
  {
#ifdef USE_OPENMP
#pragma omp for
#endif
    for (i = update_r.xmin; i < update_r.xmax; i++) {
      for (j = update_r.ymin; j < update_r.ymax; j++) {

        idx = pi(i, j, update_r.zmin);
        hz1 = &(hz_[pi(i-1, j, update_r.zmin)]);

        ey = &(ey_[idx]);
        hx = &(hx_[idx]);
        hz2 = &(hz_[idx]);

        for (k = update_r.zmin; k < update_r.zmax; k++) {
          mid = material_[idx];
          
          // Is this a plasma?
          if (vcdt_[mid] == 0.0) 
          {
            *ey = Ca_[mid] * *ey
              + Cbz_[mid] * (*hx - *(hx-1))
              + Cbx_[mid] * (*hz1 - *hz2);
          }
          else
          {
            *ey = Ca_[mid] * dy_[idx]
              + Cbz_[mid] * (*hx - *(hx-1))
              + Cbx_[mid] * (*hz1 - *hz2);
            
            dy_[idx] = *ey;
            
            *ey = dy_[idx] - sy_[idx];
            
            sy_[idx] = (1 + vcdt_[mid]) * sym1_[idx]
              - (vcdt_[mid] * sym2_[idx])
              + (omegapsq_[mid] * (1 - vcdt_[mid])) * *ey;
            
            sym2_[idx] = sym1_[idx];
            sym1_[idx] = sy_[idx];
            
            //++plasma_idx;
          }
          
          ey++;
          hx++;
          hz1++;
          hz2++;
          idx++;
        }
      }
    }
  }

}

void FreqGrid::update_ez(region_t update_r)
{
  unsigned int mid, idx, plasma_idx = 0;
  int i, j, k;
  field_t *ez, *hy1, *hy2, *hx1, *hx2;
  
  // Inner part
#ifdef USE_OPENMP
#pragma omp parallel private(mid, k, idx, ez, hy1, hy2, hx1, hx2)
#endif
  {
#ifdef USE_OPENMP
#pragma omp for
#endif
    for (i = update_r.xmin; i < update_r.xmax; i++) {
      for (j = update_r.ymin; j < update_r.ymax; j++) {
        idx = pi(i, j, update_r.zmin);
        hy2 = &(hy_[pi(i-1, j, update_r.zmin)]);
        hx1 = &(hx_[pi(i, j-1, update_r.zmin)]);

        ez = &(ez_[idx]);
        hy1 = &(hy_[idx]);
        hx2 = &(hx_[idx]);

        for (k = update_r.zmin; k < update_r.zmax; k++) {
          mid = material_[idx];
          
          // Is this a plasma?
          if (vcdt_[mid] == 0.0) 
          {
            *ez = Ca_[mid] * *ez
              + Cbx_[mid] * (*hy1 - *hy2)
              + Cby_[mid] * (*hx1 - *hx2);
          }
          else
          {
            *ez = Ca_[mid] * dz_[idx]
              + Cbx_[mid] * (*hy1 - *hy2)
              + Cby_[mid] * (*hx1 - *hx2);
            
            dz_[idx] = *ez;
            
            *ez = dz_[idx] - sz_[idx];
            
            sz_[idx] = (1 + vcdt_[mid]) * szm1_[idx]
              - (vcdt_[mid] * szm2_[idx])
              + (omegapsq_[mid] * (1 - vcdt_[mid])) * *ez;
          
            szm2_[idx] = szm1_[idx];
            szm1_[idx] = sz_[idx];
            
            //++plasma_idx;
          }
          
          ez++;
          hy1++; hy2++; hx1++; hx2++;
          idx++;
        }
      }
    }
  }
}


void FreqGrid::setup_subdomain_data(SubdomainBc *sd, Face face)
{
  Grid::setup_subdomain_data(sd, face);

  RxTxData rxtx;
  unsigned int idx_rx = 0, idx_tx = 0;
 
  rxtx.set_field_type(E);
  
  switch (face)
  {
  case FRONT: // x=dimx...
    idx_rx = pi(info_.dimx_ - 1, 0, 0);
    idx_tx = pi(info_.dimx_ - 2, 0, 0);
    break;
  case TOP: // z=dimz
    idx_rx = pi(0, 0, info_.dimz_ - 1);
    idx_tx = pi(0, 0, info_.dimz_ - 2);
    break;
  case RIGHT: // y=dimy
    idx_rx = pi(0, info_.dimy_ - 1, 0);
    idx_tx = pi(0, info_.dimy_ - 2, 0);
    break;
  case BACK: // x=0
    idx_rx = pi(0, 0, 0);
    idx_tx = pi(1, 0, 0);
    break;
  case BOTTOM: // z=0
    idx_rx = pi(0, 0, 0);
    idx_tx = pi(0, 0, 1);
    break;
  case LEFT: // y=0
    idx_rx = pi(0, 0, 0);
    idx_tx = pi(0, 1, 0);
    break;
  }

  MPI_Datatype t = get_plane_dt(face);
  rxtx.set_datatype(t);

  rxtx.set_tx_ptr(&(dx_[idx_tx]));
  rxtx.set_rx_ptr(&(dx_[idx_rx]));
  sd->add_tx_rx_data(rxtx);

  rxtx.set_tx_ptr(&(dy_[idx_tx]));
  rxtx.set_rx_ptr(&(dy_[idx_rx]));
  sd->add_tx_rx_data(rxtx);

  rxtx.set_tx_ptr(&(dz_[idx_tx]));
  rxtx.set_rx_ptr(&(dz_[idx_rx]));
  sd->add_tx_rx_data(rxtx);
  
  rxtx.set_tx_ptr(&(sx_[idx_tx]));
  rxtx.set_rx_ptr(&(sx_[idx_rx]));
  sd->add_tx_rx_data(rxtx);

  rxtx.set_tx_ptr(&(sy_[idx_tx]));
  rxtx.set_rx_ptr(&(sy_[idx_rx]));
  sd->add_tx_rx_data(rxtx);

  rxtx.set_tx_ptr(&(sz_[idx_tx]));
  rxtx.set_rx_ptr(&(sz_[idx_rx]));
  sd->add_tx_rx_data(rxtx);

  rxtx.set_tx_ptr(&(sxm1_[idx_tx]));
  rxtx.set_rx_ptr(&(sxm1_[idx_rx]));
  sd->add_tx_rx_data(rxtx);

  rxtx.set_tx_ptr(&(sym1_[idx_tx]));
  rxtx.set_rx_ptr(&(sym1_[idx_rx]));
  sd->add_tx_rx_data(rxtx);

  rxtx.set_tx_ptr(&(szm1_[idx_tx]));
  rxtx.set_rx_ptr(&(szm1_[idx_rx]));
  sd->add_tx_rx_data(rxtx);

  rxtx.set_tx_ptr(&(sxm2_[idx_tx]));
  rxtx.set_rx_ptr(&(sxm2_[idx_rx]));
  sd->add_tx_rx_data(rxtx);

  rxtx.set_tx_ptr(&(sym2_[idx_tx]));
  rxtx.set_rx_ptr(&(sym2_[idx_rx]));
  sd->add_tx_rx_data(rxtx);

  rxtx.set_tx_ptr(&(szm2_[idx_tx]));
  rxtx.set_rx_ptr(&(szm2_[idx_rx]));
  sd->add_tx_rx_data(rxtx);
  
}
