#include "FreqGrid.hh"
#include "Exceptions.hh"

#include <string.h> // for memset()
#include <iostream>

using namespace std;

FreqGrid::FreqGrid()
  : vcdt_(0), omegapsq_(0), dx_(0), sx_(0), sxm1_(0), sxm2_(0),
    dy_(0), sy_(0), sym1_(0), sym2_(0),
    dz_(0), sz_(0), szm1_(0), szm2_(0)
{}

FreqGrid::~FreqGrid()
{
  free_grid();
}

void FreqGrid::alloc_grid()
{
  cout << "Allocating grid for FreqGrid!" << endl;

  Grid::alloc_grid();

  unsigned int sz = 0;

  sz = get_ldx() * get_ldy() * get_ldz();

  // UGLY!
  unsigned int num_plasma = 0;

  // CAN'T do this because material_ hasn't been set to actual material numbers yet!
//   for (unsigned int i = 0; i < sz; i++)
//   {
//     if (vcdt_[material_[i]] != 0.0)
//       num_plasma++;
//   }

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

  int index = 0;

  vector<Material>::iterator iter = matlib.get_material_iter_begin();
  vector<Material>::iterator iter_e = matlib.get_material_iter_end();

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
    }
    ++index;
    ++iter;
  }
}

void FreqGrid::free_material()
{
  Grid::free_material();

  if (vcdt_)
    delete[] vcdt_;

  if (omegapsq_)
    delete[] omegapsq_;
}

void FreqGrid::update_ex()
{
  unsigned int mid, idx, idx2, plasma_idx = 0;
  int i, j, k;
  field_t *ex, *hz1, *hz2, *hy;

  // Inner part
  for (i = update_r_.xmin; i < update_r_.xmax; i++) {
    for (j = update_r_.ymin + 1; j < update_r_.ymax; j++) {
      
      idx = pi(i, j, update_r_.zmin + 1);
      idx2 = pi(i, j-1, update_r_.zmin + 1);
      ex = &(ex_[idx]);
      hz1 = &(hz_[idx]);
      hz2 = &(hz_[idx2]);
      hy = &(hy_[idx]);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (k = update_r_.zmin + 1; k < update_r_.zmax; k++) {
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
          *ex = Ca_[mid] * dx_[plasma_idx]
            + Cby_[mid] * (*hz1 - *hz2)
            + Cbz_[mid] * (*(hy - 1) - *hy);

          dx_[plasma_idx] = *ex;

          *ex = dx_[plasma_idx] - sx_[plasma_idx];
          
          sx_[plasma_idx] = (1 + vcdt_[mid]) * sxm1_[plasma_idx]
            - vcdt_[mid] * sxm2_[plasma_idx]
            + omegapsq_[mid] * (1 - vcdt_[mid]) * *ex;
          
          sxm2_[plasma_idx] = sxm1_[plasma_idx];
          sxm1_[plasma_idx] = sx_[plasma_idx];

          ++plasma_idx;
        }

        ex++;
        hz1++;
        hz2++;
        hy++;
      }
    }
  }
}

void FreqGrid::update_ey()
{
  unsigned int mid, idx, plasma_idx = 0;
  int i, j, k;
  field_t *ey, *hx, *hz1, *hz2;

  // Inner part
  for (i = update_r_.xmin + 1; i < update_r_.xmax; i++) {
    for (j = update_r_.ymin; j < update_r_.ymax; j++) {

      idx = pi(i, j, update_r_.zmin + 1);
      ey = &(ey_[idx]);
      hx = &(hx_[idx]);
      hz1 = &(hz_[pi(i-1, j, update_r_.zmin + 1)]);
      hz2 = &(hz_[idx]);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (k = update_r_.zmin + 1; k < update_r_.zmax; k++) {
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
          *ey = Ca_[mid] * dy_[plasma_idx]
            + Cbz_[mid] * (*hx - *(hx-1))
            + Cbx_[mid] * (*hz1 - *hz2);

          dy_[plasma_idx] = *ey;

          *ey = dy_[plasma_idx] - sy_[plasma_idx];
          
          sy_[plasma_idx] = (1 + vcdt_[mid]) * sym1_[plasma_idx]
            - vcdt_[mid] * sym2_[plasma_idx]
            + omegapsq_[mid] * (1 - vcdt_[mid]) * *ey;
          
          sym2_[plasma_idx] = sym1_[plasma_idx];
          sym1_[plasma_idx] = sy_[plasma_idx];

          ++plasma_idx;
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

void FreqGrid::update_ez()
{
  unsigned int mid, idx, plasma_idx = 0;
  int i, j, k;
  field_t *ez, *hy1, *hy2, *hx1, *hx2;
  
  // Inner part
  for (i = update_r_.xmin + 1; i < update_r_.xmax; i++) {
    for (j = update_r_.ymin + 1; j < update_r_.ymax; j++) {

      idx = pi(i, j, update_r_.zmin);
      ez = &(ez_[idx]);
      hy1 = &(hy_[idx]);
      hy2 = &(hy_[pi(i-1, j, update_r_.zmin)]);
      hx1 = &(hx_[pi(i, j-1, update_r_.zmin)]);
      hx2 = &(hx_[idx]);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (k = update_r_.zmin; k < update_r_.zmax; k++) {
        mid = material_[pi(i, j, k)];

        // Is this a plasma?
        if (vcdt_[mid] == 0.0) 
        {
          *ez = Ca_[mid] * *ez
            + Cbx_[mid] * (*hy1 - *hy2)
            + Cby_[mid] * (*hx1 - *hx2);
        }
        else
        {
          *ez = Ca_[mid] * dz_[plasma_idx]
            + Cbx_[mid] * (*hy1 - *hy2)
            + Cby_[mid] * (*hx1 - *hx2);
          
          dz_[plasma_idx] = *ez;

          *ez = dz_[plasma_idx] - sz_[plasma_idx];
          
          sz_[plasma_idx] = (1 + vcdt_[mid]) * szm1_[plasma_idx]
            - vcdt_[mid] * szm2_[plasma_idx]
            + omegapsq_[mid] * (1 - vcdt_[mid]) * *ez;
          
          szm2_[plasma_idx] = szm1_[plasma_idx];
          szm1_[plasma_idx] = sz_[plasma_idx];

          ++plasma_idx;
        }

        ez++;
        hy1++; hy2++; hx1++; hx2++;
        idx++;
      }
    }
  }
}
