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
    face_bc_[i] = EWALL;
    face_rank_[i] = 0;
  }
}

Grid::~Grid()
{
  free_grid();
  free_material();
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
  if (Ca_) {
    delete[] Ca_;
    delete[] Cb1_;
    delete[] Cb2_;
    
    delete[] Da_;
    delete[] Db1_;
    delete[] Db2_;

    Ca_ = Da_ = Cb1_ = Cb2_ = Db1_ = Db2_ = 0;
  }
}


Grid::init_datatypes()
{
  MPI_Type_contiguous(dimz_, GRID_MPI_TYPE, &z_vector_);
  MPI_Type_commit(&z_vector);

  MPI_Type_vector(dimy_, 1, dimz_, GRID_MPI_TYPE, &y_vector_);
  MPI_Type_commit(&y_vector);
  
  MPI_Type_vector(dimx_, 1, dimy_ * dimz_, GRID_MPI_TYPE, &x_vector_);
  MPI_Type_commit(&x_vector);

  // Not 100% sure about these:
  MPI_Type_vector(dimy_, 1, 0, z_vector_, &yz_plane_);
  MPI_Type_commit(&yz_plane_);

  MPI_Type_vector(dimx_, 1, dimy_, z_vector_, &xz_plane_);
  MPI_Type_commit(&xz_plane);

  MPI_Type_vector(dimx_, 1, 0, y_vector_, &xy_plane_);
  MPI_Type_commit(&xy_plane_);
}


Grid::alloc_grid()
{

  ex_ = new **field_t(dimx_);
  ey_ = new **field_t(dimx_);
  ez_ = new **field_t(dimx_);

  hx_ = new **field_t(dimx_);
  hy_ = new **field_t(dimx_);
  hz_ = new **field_t(dimx_);

  for (unsigned int i = 0; i < dimx_; i++) {
    ex_[i] = new *field_t(dimy_);
    ey_[i] = new *field_t(dimy_);
    ez_[i] = new *field_t(dimy_);

    hx_[i] = new *field_t(dimy_);
    hy_[i] = new *field_t(dimy_);
    hz_[i] = new *field_t(dimy_);

    for (unsigned int j = 0; j < dimy_; j++) {
      ex_[i][j] = new field_t(dimz_);
      ey_[i][j] = new field_t(dimz_);
      ez_[i][j] = new field_t(dimz_);

      hx_[i][j] = new field_t(dimz_);
      hy_[i][j] = new field_t(dimz_);
      hz_[i][j] = new field_t(dimz_);
    }
  }
}


Grid::load_materials(MaterialLib &matlib)
{
  
}


Grid::setup_grid(int global_x, int global_y, int global_z, 
		 int x, int y, int z)
{
  global_dimx_ = global_x; 
  global_dimy_ = global_y; 
  global_dimz_ = global_z; 

  dimx_ = x;
  dimy_ = y;
  dimz_ = z;

  alloc_grid();
}
