#include "Grid.hh"

Grid::Grid() 
  : global_dimx_(0), global_dimy_(0), global_dimz_(0), 
    dimx_(0), dimy_(0), dimz_(0), 
    Ca_(0), Cb1_(0), Cb2_(0), 
    Da_(0), Db1_(0), Db2_(0), 
    ex_(0), ey_(0), ez_(0), hx_(0), hy_(0), hz_(0), 
    ex_sum_(0), ey_sum_(0), ez_sum_(0)
{
  for (int i = 0; i < 6; i++) {
    faces_[i] = 0;
  }

  MPI_
}

Grid::~Grid()
{
  free_grid();
}
    

Grid::free_grid()
{
  // Slightly dangerous, but if one is allocated then all should 
  // be allocated. 
  if (ex_ || ey_ || ez_ || hx_ || hy_ || hz_) {
    for (unsigned int i = 0; i < dimx_; i++) {
      for (unsigned int j = 0; j < dimy_; j++) {
	delete[] ex_[i][j];
	delete[] ey_[i][j];
	delete[] ez_[i][j];

	delete[] hx_[i][j];
	delete[] hy_[i][j];
	delete[] hz_[i][j];
      }
      delete[] ex_[i];
      delete[] ey_[i];
      delete[] ez_[i];
      
      delete[] hx_[i];
      delete[] hy_[i];
      delete[] hz_[i];
    }

    delete[] ex_;
    delete[] ey_;
    delete[] ez_;
    
    delete[] hx_;
    delete[] hy_;
    delete[] hz_;
  }

  // The sums are done seperatly because they aren't always needed. 
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


Grid::free_material()
{
  
}
