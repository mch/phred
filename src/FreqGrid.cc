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
  
void FreqGrid::update_fields()
{

}

void FreqGrid::load_materials(MaterialLib &matlib)
{
  Grid::load_materials(matlib);
}
