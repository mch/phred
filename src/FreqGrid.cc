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

  ex_sum_ = new field_t **[get_ldx()];
  ey_sum_ = new field_t **[get_ldx()];
  ez_sum_ = new field_t **[get_ldx()];

  for (unsigned int i = 0; i < get_ldx(); i++) {
    ex_sum_[i] = new field_t *[get_ldy()];
    ey_sum_[i] = new field_t *[get_ldy()];
    ez_sum_[i] = new field_t *[get_ldy()];

    for (unsigned int j = 0; j < get_ldy(); j++) {  
      ex_sum_[i][j] = new field_t[get_ldz()];
      ey_sum_[i][j] = new field_t[get_ldz()];
      ez_sum_[i][j] = new field_t[get_ldz()];

    }
  }
}

void FreqGrid::free_grid()
{
  Grid::free_grid();

  if (ex_sum_ || ey_sum_ || ez_sum_) {
    for (unsigned int i = 0; i < get_ldx(); i++) {
      for (unsigned int j = 0; j < get_ldy(); j++) {  
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
