#include "Grid.hh"

Grid::Grid() 
  : global_dimx_(0), global_dimy_(0), global_dimz_(0), 
    dimx_(0), dimy_(0), dimz_(0), 
    deltax_(0), deltay_(0), deltaz_(0), deltat_(0), 
    num_materials_(0),
    Ca_(0), Cbx_(0), Cby_(0), Cbz_(0),
    Da_(0), Dbx_(0), Dby_(0), Dbz_(0),
    ex_(0), ey_(0), ez_(0), hx_(0), hy_(0), hz_(0), 
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
    

void Grid::free_grid()
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
}


void Grid::free_material()
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


void Grid::init_datatypes()
{
  MPI_Type_contiguous(dimz_, GRID_MPI_TYPE, &z_vector_);
  MPI_Type_commit(&z_vector_);

  MPI_Type_vector(dimy_, 1, dimz_, GRID_MPI_TYPE, &y_vector_);
  MPI_Type_commit(&y_vector_);
  
  MPI_Type_vector(dimx_, 1, dimy_ * dimz_, GRID_MPI_TYPE, &x_vector_);
  MPI_Type_commit(&x_vector_);

  // Not 100% sure about these:
  MPI_Type_vector(dimy_, 1, 0, z_vector_, &yz_plane_);
  MPI_Type_commit(&yz_plane_);

  MPI_Type_vector(dimx_, 1, dimy_, z_vector_, &xz_plane_);
  MPI_Type_commit(&xz_plane_);

  MPI_Type_vector(dimx_, 1, 0, y_vector_, &xy_plane_);
  MPI_Type_commit(&xy_plane_);
}


void Grid::alloc_grid()
{

  ex_ = new field_t **[dimx_];
  ey_ = new field_t **[dimx_];
  ez_ = new field_t **[dimx_];

  hx_ = new field_t **[dimx_];
  hy_ = new field_t **[dimx_];
  hz_ = new field_t **[dimx_];

  material_ = new unsigned int **[dimx_];

  for (unsigned int i = 0; i < dimx_; i++) {
    ex_[i] = new field_t *[dimy_];
    ey_[i] = new field_t *[dimy_];
    ez_[i] = new field_t *[dimy_];

    hx_[i] = new field_t *[dimy_];
    hy_[i] = new field_t *[dimy_];
    hz_[i] = new field_t *[dimy_];

    material_[i] = new unsigned int *[dimy_];

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


void Grid::load_materials(MaterialLib &matlib)
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


void Grid::setup_grid(unsigned int global_x, unsigned int global_y, 
                      unsigned int global_z, 
                      unsigned int start_x, unsigned int start_y, 
                      unsigned int start_z, 
                      unsigned int x, unsigned int y, unsigned int z, 
                      delta_t deltax, delta_t deltay, delta_t deltaz,
                      delta_t deltat)
{
  global_dimx_ = global_x; 
  global_dimy_ = global_y; 
  global_dimz_ = global_z; 

  start_x_ = start_x;
  start_y_ = start_y;
  start_z_ = start_z;

  dimx_ = x;
  dimy_ = y;
  dimz_ = z;

  deltax_ = deltax;
  deltay_ = deltay;
  deltaz_ = deltaz;
  deltat_ = deltat;

  alloc_grid();
}


void Grid::define_box(unsigned int x_start, unsigned int x_stop, 
                      unsigned int y_start, unsigned int y_stop, 
                      unsigned int z_start, unsigned int z_stop, 
                      unsigned int mat_index)
{
  // Given coordinates are global, so we have to convert them to local. 
  unsigned int xs, ys, zs, xe, ye, ze;

  xs = (start_x_ > x_start) ? start_x_ : x_start - start_x_;
  ys = (start_y_ > y_start) ? start_y_ : y_start - start_y_;
  zs = (start_z_ > z_start) ? start_z_ : z_start - start_z_;

  xe = (start_x_ + dimx_ > x_stop) ? start_x_ + dimx_ : x_stop - start_x_;
  ye = (start_y_ + dimy_ > y_stop) ? start_y_ + dimy_ : y_stop - start_y_;
  ze = (start_z_ + dimz_ > z_stop) ? start_z_ + dimz_ : z_stop - start_z_;

  for (unsigned int i = xs; i < xe; i++)
  {
    for (unsigned int j = ys; j < ye; j++)
    {
      for (unsigned int k = zs; k < ze; k++)
      {
        material_[i][j][k] = mat_index;
      }
    }
  }
}

void Grid::set_boundary(unsigned int face, BoundaryCondition bc)
{
  face_bc_[face] = bc;
}

void Grid::set_face_rank(unsigned int face, int rank)
{
  face_rank_[face] = rank;
}


// Straight out of Taflove.
void Grid::update_fields()
{
  unsigned int mid;
  for (unsigned int i = 1; i < dimx_ - 1; i++) {
    for (unsigned int j = 1; j < dimy_ - 1; j++) {
      for (unsigned int k = 1; k < dimz_ - 1; k++) {
        mid = material_[i][j][k];

        // Electric
        ex_[i][j][k] = Ca_[mid] * ex_[i][j][k]
          + Cby_[mid] * (hz_[i][j][k] - hz_[i][j-1][k])
          + Cbz_[mid] * (hy_[i][j][k-1] - hy_[i][j][k]);

        ey_[i][j][k] = Ca_[mid] * ey_[i][j][k]
          + Cbz_[mid] * (hx_[i][j][k] - hx_[i][j][k-1])
          + Cbx_[mid] * (hz_[i-1][j][k] - hz_[i][j][k]);

        ez_[i][j][k] = Ca_[mid] * ez_[i][j][k]
          + Cbx_[mid] * (hy_[i][j][k] - hy_[i-1][j][k])
          + Cby_[mid] * (hx_[i][j-1][k] - hx_[i][j][k]);          

        // Magnetic
        hx_[i][j][k] = Da_[mid] * hx_[i][j][k]
          + Dby_[mid] * (ez_[i][j][k] - ez_[i][j+1][k])
          + Dbz_[mid] * (ey_[i][j][k+1] - ey_[i][j][k]);

        hy_[i][j][k] = Da_[mid] * hy_[i][j][k]
          + Dbz_[mid] * (ex_[i][j][k] - ex_[i][j][k+1])
          + Dbx_[mid] * (ez_[i+1][j][k] - ez_[i][j][k]);

        hz_[i][j][k] = Da_[mid] * hz_[i][j][k]
          + Dbx_[mid] * (ey_[i][j][k] - ey_[i+1][j][k])
          + Dby_[mid] * (ex_[i][j+1][k] - ex_[i][j][k]);
      }
    }
  }
}

void Grid::apply_boundaries()
{

}

// Most of this should be inline, in the header file!

unsigned int Grid::get_gdx()
{
  return global_dimx_;
}

unsigned int Grid::get_gdy()
{
  return global_dimy_;
}

unsigned int Grid::get_gdz()
{
  return global_dimz_;
}


unsigned int Grid::get_lsx()
{
  return start_x_;
}

unsigned int Grid::get_lsy()
{
  return start_y_;
}

unsigned int Grid::get_lsz()
{
  return start_x_;
}

unsigned int Grid::get_ldx()
{
  return dimx_;
}

unsigned int Grid::get_ldy()
{
  return dimy_;
}

unsigned int Grid::get_ldz()
{
  return dimz_;
}

delta_t Grid::get_deltax()
{
  return deltax_;
}

delta_t Grid::get_deltay()
{
  return deltay_;
}

delta_t Grid::get_deltaz()
{
  return deltaz_;
}

delta_t Grid::get_deltat()
{
  return deltat_;
}

void Grid::set_ex(unsigned int x, unsigned int y, 
                  unsigned int z, field_t val)
{
  ex_[x][y][z] = val;
}

void Grid::set_ey(unsigned int x, unsigned int y, 
                  unsigned int z, field_t val)
{
  ey_[x][y][z] = val;
}

void Grid::set_ez(unsigned int x, unsigned int y, 
                  unsigned int z, field_t val)
{
  ez_[x][y][z] = val;
}

void Grid::set_hx(unsigned int x, unsigned int y, 
                  unsigned int z, field_t val)
{
  hx_[x][y][z] = val;
}

void Grid::set_hy(unsigned int x, unsigned int y, 
                  unsigned int z, field_t val)
{
  hy_[x][y][z] = val;
}

void Grid::set_hz(unsigned int x, unsigned int y, 
                  unsigned int z, field_t val)
{
  hz_[x][y][z] = val;
}

