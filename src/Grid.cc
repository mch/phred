#include "Grid.hh"

Grid::Grid() 
  : global_dimx_(0), global_dimy_(0), global_dimz_(0), 
    dimx_(0), dimy_(0), dimz_(0), 
    deltax_, deltay_, deltaz_, deltat_, 
    num_materials_(0),
    Ca_(0), Cbx_(0), Cby_(0), Cbz_(0),
    Da_(0), Dbx_(0), Dby_(0), Dbz_(0),
    ex_(0), ey_(0), ez_(0), hx_(0), hy_(0), hz_(0), 
    ex_sum_(0), ey_sum_(0), ez_sum_(0),
    material_(0)
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
    delete[] Cbx_;
    delete[] Cby_;
    delete[] Cbz_;
    
    delete[] Da_;
    delete[] Dbx_;
    delete[] Dby_;
    delete[] Dbz_;

    Ca_ = Da_ = Cbx_ = Cby_ = Cbz_ = Dbx_ = Dby_ = Dbz_ = 0;
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

  material_ = new **unsigned int(dimx_);

  for (unsigned int i = 0; i < dimx_; i++) {
    ex_[i] = new *field_t(dimy_);
    ey_[i] = new *field_t(dimy_);
    ez_[i] = new *field_t(dimy_);

    hx_[i] = new *field_t(dimy_);
    hy_[i] = new *field_t(dimy_);
    hz_[i] = new *field_t(dimy_);

    material_[i] = new *unsigned int(dimy_);

    for (unsigned int j = 0; j < dimy_; j++) {
      ex_[i][j] = new field_t(dimz_);
      ey_[i][j] = new field_t(dimz_);
      ez_[i][j] = new field_t(dimz_);

      hx_[i][j] = new field_t(dimz_);
      hy_[i][j] = new field_t(dimz_);
      hz_[i][j] = new field_t(dimz_);

      material_[i][j] = new unsigned int(dimz_);

      for (unsigned int k = 0; k < dimz_; k++)
      {
        ex_[i][j][k] = 0;
        ey_[i][j][k] = 0;
        ez_[i][j][k] = 0;

        hx_[i][j][k] = 0;
        hy_[i][j][k] = 0;
        hz_[i][j][k] = 0;

        material_[i][j][k] = 0; // PEC
      }
    }
  }
}


Grid::load_materials(MaterialLib &matlib)
{
  // Clear up any material data that may already be loaded
  free_material();

  int num_mat = matlib.num_materials();
  Ca_ = new mat_coef_t(num_mat);
  Da_ = new mat_coef_t(num_mat);

  Cbx_ = new mat_coef_t(num_mat);
  Dbx_ = new mat_coef_t(num_mat);

  // Save some memory if possible. 
  if (deltay_ == deltax_)
  {
    Cby_ = Cbx_;
    Dby_ = Dbx_;
  } else {
    Cby_ = new mat_coef_t(num_mat);
    Dby_ = new mat_coef_t(num_mat);
  }

  if (deltaz_ == deltax_)
  {
    Cbz_ = Cbx_;
    Dbz_ = Dbx_;
  } else if (deltaz_ == deltay_) {
    Cbz_ = Cby_;
    Dbz_ = Dby_;
  } else {
    Cbz_ = new mat_coef_t(num_mat);
    Dbz_ = new mat_coef_t(num_mat);
  }

  vector<Material>::iterator iter = matlib.get_material_iter_begin();
  vector<Material>::iterator iter_e = matlib.get_material_iter_end();

  // The first one is always PEC
  int index = 0;

  Ca_[index] = 1;
  Cbx_[index] = Cby_[index] = Cbz_[index] = 0;

  Da_[index] = 1;
  Dbx_[index] = Dby_[index] = Dbz_[index] = 0;
  
  ++index;
  while (iter != iter_e) 
  {
    // Make the code cleaner with short var names
    mat_prop_t eps = (*iter).get_epsilon();
    mat_prop_t sig = (*iter).get_sigma();
    mat_prop_t mu = (*iter).get_mu();
    mat_prop_t sigs = (*iter).get_sigma_star();

    if (eps == 0 || mu == 0)
    {
      cerr << "Something is wrong with the material library:\n" 
           << " -> Material cannot have permittivities or permeabilities\n"
           << "    of zero. Perfect electric conductor can have eps=0, \n"
           << "    but that is a special material defined by phred.\n\n"
           << "Program aborting. Check material library." << endl;
      exit(1);
    }
    
    Ca_[index] = (1 - (sig * deltat_ * 0.5)/eps) / 
                 (1 + (sig * deltat_ * 0.5)/eps);

    Da_[index] = (1 - (sigs * deltat_ * 0.5)/mu) / 
                 (1 + (sigs * deltat_ * 0.5)/mu);

    
    Cbx_[index] = (deltat_ / (eps * deltax_)) / 
                  (1 + (sig * deltat_ * 0.5)/eps);

    Dbx_[index] = (deltat_ / (mu * deltax_)) / 
                  (1 + (sigs * deltat_ * 0.5)/mu);

    if (deltay_ != deltax_)
    {    
      Cby_[index] = (deltat_ / (eps * deltay_)) / 
                    (1 + (sig * deltat_ * 0.5)/eps);

      Dby_[index] = (deltat_ / (mu * deltay_)) / 
                    (1 + (sigs * deltat_ * 0.5)/mu);
    }

    if (deltaz_ != deltax_ && deltaz_ != deltay_)
    {
      Cbz_[index] = (deltat_ / (eps * deltaz_)) / 
                    (1 + (sig * deltat_ * 0.5)/eps);

      Dbz_[index] = (deltat_ / (mu * deltaz_)) / 
                    (1 + (sigs * deltat_ * 0.5)/mu);
    }
    
  }
}


void Grid::setup_grid(int global_x, int global_y, int global_z, 
                      int x, int y, int z, 
                      delta_t deltax, delta_t deltay, delta_t deltaz,
                      delta_t deltat)
{
  global_dimx_ = global_x; 
  global_dimy_ = global_y; 
  global_dimz_ = global_z; 

  dimx_ = x;
  dimy_ = y;
  dimz_ = z;

  deltax_ = deltax;
  deltay_ = deltay;
  deltaz_ = deltaz;
  deltat_ = deltat;

  alloc_grid();
}


void Grid::define_box(int x_start, int x_stop, int y_start, int y_stop, 
                      int z_start, int z_stop, unsigned int mat_index)
{
  // Given coordinates are global, so we have to convert them to local. 
}

void Grid::set_boundary(unsigned int face, BoundaryCondition bc)
{
  face_bc_[face] = bc;
}

void Grid::set_face_rank(unsigned int face, int rank)
{
  face_rank_[face] = rank;
}
