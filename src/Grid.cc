#include "Grid.hh"

#include "Ewall.hh"
#include "Hwall.hh"

Grid::Grid() 
  : num_materials_(0),
    Ca_(0), Cbx_(0), Cby_(0), Cbz_(0),
    Da_(0), Dbx_(0), Dby_(0), Dbz_(0),
    ex_(0), ey_(0), ez_(0), hx_(0), hy_(0), hz_(0), 
    material_(0), define_(true)
{

}

Grid::~Grid()
{
  define_ = true;

  free_grid();
  free_material();
}
    
void Grid::set_define_mode(bool d)
{
  bool ok = true;

  if (!d) 
  {
    // Sanity checks
    
    // Stability check
    
    // Calculate update region_t by considering the thickness of the PML's. 
    update_r_.xmin = 0;
    update_r_.xmax = info_.dimx_;
    update_r_.ymin = 0;
    update_r_.ymax = info_.dimy_;
    update_r_.zmin = 0;
    update_r_.zmax = info_.dimz_;

    // The domain decomp algorithm will take care of assigning the
    // boundary conditions sensibly, so we don't have to worry about
    // wether or not we are really on a boundary that has thickness
    // (i.e. a PML face)
    unsigned int thickness = 0;
    for (int i = 0; i < 6; i++)
    {
      thickness = info_.get_face_thickness(static_cast<Face>(i));

      switch (i) {
      case FRONT:
        update_r_.xmax -= thickness;
      case BACK:
        update_r_.xmin += thickness;
      case TOP:
        update_r_.zmax -= thickness;
      case BOTTOM:
        update_r_.zmin += thickness;
      case LEFT:
        update_r_.ymin += thickness;
      case RIGHT:
        update_r_.ymax -= thickness;
      }
    }

    if (ok)
      define_ = d;
    else 
    {
      cerr << "The grid is not in a sane condition which can be reasonably be solved. " << endl;
    }
  } else {
    define_ = d;
  }
}

void Grid::free_grid()
{
  if (!define_)
  {
    cerr << "Unable to free grid data; the grid is not in define mode." << endl;
    return;
  }

  // Slightly dangerous, but if one is allocated then all should 
  // be allocated. 
  if (ex_ || ey_ || ez_ || hx_ || hy_ || hz_) {
    for (unsigned int i = 0; i < get_ldx(); i++) {
      for (unsigned int j = 0; j < get_ldy(); j++) {
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
  if (!define_)
  {
    cerr << "Unable to free material data; the grid is not in define mode." << endl;
    return;
  }

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
  MPI_Type_contiguous(get_ldz(), GRID_MPI_TYPE, &z_vector_);
  MPI_Type_commit(&z_vector_);

  MPI_Type_vector(get_ldy(), 1, get_ldz(), GRID_MPI_TYPE, &y_vector_);
  MPI_Type_commit(&y_vector_);
  
  MPI_Type_vector(get_ldx(), 1, get_ldy() * get_ldz(), GRID_MPI_TYPE, &x_vector_);
  MPI_Type_commit(&x_vector_);

  // Not 100% sure about these:
  MPI_Type_vector(get_ldy(), 1, 0, z_vector_, &yz_plane_);
  MPI_Type_commit(&yz_plane_);

  MPI_Type_vector(get_ldx(), 1, get_ldy(), z_vector_, &xz_plane_);
  MPI_Type_commit(&xz_plane_);

  MPI_Type_vector(get_ldx(), 1, 0, y_vector_, &xy_plane_);
  MPI_Type_commit(&xy_plane_);
}


void Grid::alloc_grid()
{
  if (!define_)
  {
    cerr << "Unable to allocate grid data; the grid is not in define mode." << endl;
    return;
  }

  ex_ = new field_t **[get_ldx()];
  ey_ = new field_t **[get_ldx()];
  ez_ = new field_t **[get_ldx()];

  hx_ = new field_t **[get_ldx()];
  hy_ = new field_t **[get_ldx()];
  hz_ = new field_t **[get_ldx()];

  material_ = new unsigned int **[get_ldx()];

  for (unsigned int i = 0; i < get_ldx(); i++) {
    ex_[i] = new field_t *[get_ldy()];
    ey_[i] = new field_t *[get_ldy()];
    ez_[i] = new field_t *[get_ldy()];

    hx_[i] = new field_t *[get_ldy()];
    hy_[i] = new field_t *[get_ldy()];
    hz_[i] = new field_t *[get_ldy()];

    material_[i] = new unsigned int *[get_ldy()];

    for (unsigned int j = 0; j < get_ldy(); j++) {
      ex_[i][j] = new field_t[get_ldz()];
      ey_[i][j] = new field_t[get_ldz()];
      ez_[i][j] = new field_t[get_ldz()];

      hx_[i][j] = new field_t[get_ldz()];
      hy_[i][j] = new field_t[get_ldz()];
      hz_[i][j] = new field_t[get_ldz()];

      material_[i][j] = new unsigned int[get_ldz()];

      for (unsigned int k = 0; k < get_ldz(); k++)
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
  if (!define_)
  {
    cerr << "Unable to load material data; the grid is not in define mode." << endl;
    return;
  }

  // Clear up any material data that may already be loaded
  free_material();

  int num_mat = matlib.num_materials();
  Ca_ = new mat_coef_t(num_mat);
  Da_ = new mat_coef_t(num_mat);

  Cbx_ = new mat_coef_t(num_mat);
  Dbx_ = new mat_coef_t(num_mat);

  // Save some memory if possible. 
  if (get_deltay() == get_deltax())
  {
    Cby_ = Cbx_;
    Dby_ = Dbx_;
  } else {
    Cby_ = new mat_coef_t(num_mat);
    Dby_ = new mat_coef_t(num_mat);
  }

  if (get_deltaz() == get_deltax())
  {
    Cbz_ = Cbx_;
    Dbz_ = Dbx_;
  } else if (get_deltaz() == get_deltay()) {
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
    
    Ca_[index] = (1 - (sig * get_deltat() * 0.5)/eps) / 
                 (1 + (sig * get_deltat() * 0.5)/eps);

    Da_[index] = (1 - (sigs * get_deltat() * 0.5)/mu) / 
                 (1 + (sigs * get_deltat() * 0.5)/mu);

    
    Cbx_[index] = (get_deltat() / (eps * get_deltax())) / 
                  (1 + (sig * get_deltat() * 0.5)/eps);

    Dbx_[index] = (get_deltat() / (mu * get_deltax())) / 
                  (1 + (sigs * get_deltat() * 0.5)/mu);

    if (get_deltay() != get_deltax())
    {    
      Cby_[index] = (get_deltat() / (eps * get_deltay())) / 
                    (1 + (sig * get_deltat() * 0.5)/eps);

      Dby_[index] = (get_deltat() / (mu * get_deltay())) / 
                    (1 + (sigs * get_deltat() * 0.5)/mu);
    }

    if (get_deltaz() != get_deltax() && get_deltaz() != get_deltay())
    {
      Cbz_[index] = (get_deltat() / (eps * get_deltaz())) / 
                    (1 + (sig * get_deltat() * 0.5)/eps);

      Dbz_[index] = (get_deltat() / (mu * get_deltaz())) / 
                    (1 + (sigs * get_deltat() * 0.5)/mu);
    }

    ++iter;
    ++index;
  }
}


void Grid::setup_grid(const GridInfo &info)
{
  if (!define_)
  {
    cerr << "Unable to setup grid; the grid is not in define mode." << endl;
    return;
  }

  info_ = info;
}


void Grid::define_box(unsigned int x_start, unsigned int x_stop, 
                      unsigned int y_start, unsigned int y_stop, 
                      unsigned int z_start, unsigned int z_stop, 
                      unsigned int mat_index)
{
  // Given coordinates are global, so we have to convert them to local. 
  unsigned int xs, ys, zs, xe, ye, ze;

  if (!define_)
  {
    cerr << "Unable to define a box; the grid is not in define mode." << endl;
    return;
  }

  xs = (get_lsx() > x_start) ? get_lsx() - 1 : x_start - get_lsx();
  ys = (get_lsy() > y_start) ? get_lsy() - 1 : y_start - get_lsy();
  zs = (get_lsz() > z_start) ? get_lsz() - 1 : z_start - get_lsz();

  xe = (x_stop >= get_lsx() + get_ldx()) 
    ? get_ldx() : x_stop - get_lsx() + 1;

  ye = (y_stop >= get_lsy() + get_ldy()) 
    ? get_ldy() : y_stop - get_lsy() + 1;

  ze = (z_stop >= get_lsz() + get_ldz()) 
    ? get_ldz() : z_stop - get_lsz() + 1;

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

// Straight out of Taflove.
void Grid::update_fields()
{
  if (define_)
  {
    cerr << "Unable to update fields; the grid is in define mode." << endl;
    return;
  }

  update_ex();
  update_ey();
  update_ez();

  update_hx();
  update_hy();
  update_hz();
}

// Straight out of Taflove.
void Grid::update_ex() 
{
  unsigned int mid, i, j, k;
  
  // Inner part
  for (i = 0; i < get_ldx(); i++) {
    for (j = 1; j < get_ldy(); j++) {
      for (k = 1; k < get_ldz(); k++) {
        mid = material_[i][j][k];

        ex_[i][j][k] = Ca_[mid] * ex_[i][j][k]
          + Cby_[mid] * (hz_[i][j][k] - hz_[i][j-1][k])
          + Cbz_[mid] * (hy_[i][j][k-1] - hy_[i][j][k]);
      }
    }
  }
}

// Straight out of Taflove.
void Grid::update_ey() 
{
  unsigned int mid, i, j, k;
  
  // Inner part
  for (i = 1; i < get_ldx(); i++) {
    for (j = 0; j < get_ldy(); j++) {
      for (k = 1; k < get_ldz(); k++) {
        mid = material_[i][j][k];

        ey_[i][j][k] = Ca_[mid] * ey_[i][j][k]
          + Cbz_[mid] * (hx_[i][j][k] - hx_[i][j][k-1])
          + Cbx_[mid] * (hz_[i-1][j][k] - hz_[i][j][k]);
      }
    }
  }
}

// Straight out of Taflove.
void Grid::update_ez() 
{
  unsigned int mid, i, j, k;
  
  // Inner part
  for (i = 1; i < get_ldx(); i++) {
    for (j = 1; j < get_ldy(); j++) {
      for (k = 0; k < get_ldz(); k++) {
        mid = material_[i][j][k];

        ez_[i][j][k] = Ca_[mid] * ez_[i][j][k]
          + Cbx_[mid] * (hy_[i][j][k] - hy_[i-1][j][k])
          + Cby_[mid] * (hx_[i][j-1][k] - hx_[i][j][k]);
      }
    }
  }
}

// Straight out of Taflove.
void Grid::update_hx()
{
  unsigned int mid, i, j, k;

  for (i = 0; i < get_ldx(); i++) {
    for (j = 0; j < get_ldy() - 1; j++) {
      for (k = 0; k < get_ldz() - 1; k++) {
        mid = material_[i][j][k];

        hx_[i][j][k] = Da_[mid] * hx_[i][j][k]
          + Dby_[mid] * (ez_[i][j][k] - ez_[i][j+1][k])
          + Dbz_[mid] * (ey_[i][j][k+1] - ey_[i][j][k]);        
      }
    }
  }
}

// Straight out of Taflove.
void Grid::update_hy()
{
  unsigned int mid, i, j, k;

  for (i = 0; i < get_ldx() - 1; i++) {
    for (j = 0; j < get_ldy(); j++) {
      for (k = 0; k < get_ldz() - 1; k++) {
        mid = material_[i][j][k];

        hy_[i][j][k] = Da_[mid] * hy_[i][j][k]
          + Dbz_[mid] * (ex_[i][j][k] - ex_[i][j][k+1])
          + Dbx_[mid] * (ez_[i+1][j][k] - ez_[i][j][k]);        
      }
    }
  }
}

// Straight out of Taflove.
void Grid::update_hz()
{
  unsigned int mid, i, j, k;

  for (i = 0; i < get_ldx() - 1; i++) {
    for (j = 0; j < get_ldy() - 1; j++) {
      for (k = 0; k < get_ldz(); k++) {
        mid = material_[i][j][k];

        hz_[i][j][k] = Da_[mid] * hz_[i][j][k]
          + Dbx_[mid] * (ey_[i][j][k] - ey_[i+1][j][k])
          + Dby_[mid] * (ex_[i][j+1][k] - ex_[i][j][k]);        
      }
    }
  }
}

void Grid::apply_boundaries()
{
  if (define_)
  {
    cerr << "Unable to apply boundary conditions; the grid is in define mode." << endl;
    return;
  }

  info_.apply_boundaries(*this);
}


region_t Grid::global_to_local(region_t in)
{
  region_t r;

  r.xmin = (get_lsx() > in.xmin) ? get_lsx() - 1
    : in.xmin - get_lsx();
  r.ymin = (get_lsy() > in.ymin) ? get_lsy() - 1
    : in.ymin - get_lsy();
  r.zmin = (get_lsz() > in.zmin) ? get_lsz() - 1
    : in.zmin - get_lsz();

  r.xmax = (in.xmax >= get_lsx() + get_ldx()) 
    ? get_ldx() : in.xmax - get_lsx() + 1;

  r.ymax = (in.ymax >= get_lsy() + get_ldy()) 
    ? get_ldy() : in.ymax - get_lsy() + 1;

  r.zmax = (in.zmax >= get_lsz() + get_ldz()) 
    ? get_ldz() : in.zmax - get_lsz() + 1;

  return r;
}

region_t Grid::global_to_local(unsigned int x_start, unsigned int x_stop, 
                             unsigned int y_start, unsigned int y_stop, 
                             unsigned int z_start, unsigned int z_stop)
{
  region_t result;
  
  result.xmin = x_start;
  result.xmax = x_stop;
  result.ymin = y_start;
  result.ymax = y_stop;
  result.zmin = z_start;
  result.zmax = z_stop;

  return global_to_local(result);
}

