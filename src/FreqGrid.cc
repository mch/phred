#include "FreqGrid.hh"

FreqGrid::FreqGrid()
  : ex_sum_(0), ey_sum_(0), ez_sum_(0)
{}

FreqGrid::~FreqGrid()
{
  free_grid();
}

void FreqGrid::alloc_grid()
{
  Grid::alloc_grid();

  unsigned int sz = 0;

  sz = get_ldx() * get_ldy() * get_ldz();

  if (sz > 0) 
  {  
    ex_sum_ = new field_t[sz];
    ey_sum_ = new field_t[sz];
    ez_sum_ = new field_t[sz];
  }
}

void FreqGrid::free_grid()
{
  Grid::free_grid();

  if (ex_sum_ || ey_sum_ || ez_sum_) {
    delete[] ex_sum_;
    delete[] ey_sum_;
    delete[] ez_sum_;
  }
}
  
void FreqGrid::load_materials(MaterialLib &matlib)
{
  Grid::load_materials(matlib);
}

void FreqGrid::free_material()
{

}

void FreqGrid::update_ex()
{
 unsigned int mid, idx, idx2;
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
        
        *ex = Ca_[mid] * *ex
          + Cby_[mid] * (*hz1 - *hz2)
          + Cbz_[mid] * (*(hy - 1) - *hy);

        ex++;
        hz1++;
        hz2++;
        hy++;
      }
    }
  }
}

void FreqGrid::update_ey()
{}

void FreqGrid::update_ez()
{}

void FreqGrid::update_hx()
{}

void FreqGrid::update_hy()
{}

void FreqGrid::update_hz()
{}
