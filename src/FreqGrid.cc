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

  ex_sum_ = new field_t **[dimx_];
  ey_sum_ = new field_t **[dimx_];
  ez_sum_ = new field_t **[dimx_];

  for (unsigned int i = 0; i < dimx_; i++) {
    ex_sum_[i] = new field_t *[dimy_];
    ey_sum_[i] = new field_t *[dimy_];
    ez_sum_[i] = new field_t *[dimy_];

    for (unsigned int j = 0; j < dimy_; j++) {  
      ex_sum_[i][j] = new field_t[dimz_];
      ey_sum_[i][j] = new field_t[dimz_];
      ez_sum_[i][j] = new field_t[dimz_];

    }
  }
}

void FreqGrid::free_grid()
{
  Grid::free_grid();

  if (ex_sum_ || ey_sum_ || ez_sum_) {
    for (unsigned int i = 0; i < dimx_; i++) {
      for (unsigned int j = 0; j < dimy_; j++) {  
	delete[] ex_sum_[i][j];
	delete[] ey_sum_[i][j];
	delete[] ez_sum_[i][j];
      }
      delete[] ex_sum_[i];
      delete[] ey_sum_[i];
      delete[] ez_sum_[i];
    }
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
