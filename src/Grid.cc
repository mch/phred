#include "config.h"

#include "Grid.hh"

#include "Ewall.hh"
#include "Hwall.hh"
#include "Constants.hh"
#include "Exceptions.hh"

#include <mpi.h>

Grid::Grid() 
  : num_materials_(0),
    Ca_(0), Cbx_(0), Cby_(0), Cbz_(0),
    Da_(0), Dbx_(0), Dby_(0), Dbz_(0),
    ex_(0), ey_(0), ez_(0), hx_(0), hy_(0), hz_(0), 
    material_(0), types_alloced_(false), define_(true)
{

}

Grid::Grid(const Grid &rhs)
{
  *this = rhs;
}

const Grid &Grid::operator=(const Grid &rhs)
{
  info_ = rhs.info_;
  update_r_ = rhs.update_r_;
  num_materials_ = rhs.num_materials_;
  Ca_ = Cbx_ = Cby_ = Cbz_ = Da_ = Dbx_ = Dby_ = Dbz_ = 0;
  ex_ = ey_ = ez_ = hx_ = hy_ = hz_ = 0;
  material_ = 0;
  xy_plane_ = rhs.xy_plane_;
  yz_plane_ = rhs.yz_plane_;
  xz_plane_ = rhs.xz_plane_;
  x_vector_ = rhs.x_vector_;
  y_vector_ = rhs.y_vector_;
  z_vector_ = rhs.z_vector_;
  define_ = rhs.define_;

  return *this;
}

Grid::~Grid()
{
  define_ = true;

  free_grid();
  free_material();
  free_datatypes();
}
    
void Grid::set_define_mode(bool d)
{
  bool ok = true;

  if (!d) 
  {
    // Sanity checks (like the grid has non zero size in all dimensions)
    
    // Stability check
    float temp = C * sqrt( 1/(pow(get_deltax(), static_cast<float>(2.0))) + 
                           1/(pow(get_deltay(), static_cast<float>(2.0))) + 
                           1/(pow(get_deltaz(), static_cast<float>(2.0))));

    if (get_deltat() > 1/temp)
      throw StabilityException();

    // Initialize the MPI Derived data types
    init_datatypes();

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
        break;
      case BACK:
        update_r_.xmin += thickness;
        break;
      case TOP:
        update_r_.zmax -= thickness;
        break;
      case BOTTOM:
        update_r_.zmin += thickness;
        break;
      case LEFT:
        update_r_.ymin += thickness;
        break;
      case RIGHT:
        update_r_.ymax -= thickness;
        break;
      }

      // Initialize the PML's
      Pml *p = dynamic_cast<Pml *>(&info_.get_boundary(static_cast<Face>(i)));
      if (p) {
        p->setup(static_cast<Face>(i), *this);

        // Check ajacent faces and see if there are any subdomains
        // that need to have data shared across them.
        for (int j = 0; j < 6; j++)
        {
          if (j != i 
              && (j % 2 && (j-1) != i) // j is odd
              && ( !(j % 2) && (j+1) != i)) // j is even
          {
            if (info_.get_bc_type(static_cast<Face>(j)) == SUBDOMAIN)
            {
              SubdomainBc *sd = dynamic_cast<SubdomainBc *>
                (&info_.get_boundary(static_cast<Face>(j)));
              
              p->add_sd_bcs(sd, static_cast<Face>(i), 
                            static_cast<Face>(j));
            }
          }
        }
      }
    }
    
//     cout << "Grid update region: x={" << update_r_.xmin << "," 
//          << update_r_.xmax << "}, y={" << update_r_.ymin << "," 
//          << update_r_.ymax << "}, z={" << update_r_.zmin << ","
//          << update_r_.zmax << "}.\n";

    // Calculate common PML coefficients. 
    pml_common_.init_coeffs(*this);

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

  if (ex_)
    delete[] ex_;

  if (ey_)
    delete[] ey_;

  if (ez_)
    delete[] ez_;
    
  if (hx_)
    delete[] hx_;

  if (hy_)
    delete[] hy_;

  if (hz_)
    delete[] hz_;

  ex_ = ey_ = ez_ = hx_ = hy_ = hz_ = 0;

}


void Grid::free_material()
{
  if (!define_)
  {
    cerr << "Unable to free material data; the grid is not in define mode." << endl;
    return;
  }

  if (Ca_) 
    delete[] Ca_;

  if (Da_)
    delete[] Da_;

  if (Cbx_)
    delete[] Cbx_;

  if (Dbx_)
    delete[] Dbx_;

  if (get_deltay() != get_deltax())
  {
    if (Cby_)
      delete[] Cby_;
    
    if (Dby_)
      delete[] Dby_;
  }
  
  if (get_deltaz() != get_deltay() && get_deltaz() != get_deltax())
  {
    if (Cbz_)
      delete[] Cbz_;
    
    if (Dbz_)
      delete[] Dbz_;
  }
  
  Ca_ = Da_ = Cbx_ = Cby_ = Cbz_ = Dbx_ = Dby_ = Dbz_ = 0;
}

void Grid::free_datatypes()
{
  if (types_alloced_) {
    MPI_Type_free(&xy_plane_);
    MPI_Type_free(&xz_plane_);
    MPI_Type_free(&yz_plane_);
    
    MPI_Type_free(&x_vector_);
    MPI_Type_free(&y_vector_);
    MPI_Type_free(&z_vector_);
    
    types_alloced_ = false;
  }
}

void Grid::init_datatypes()
{
  free_datatypes();

  MPI_Type_contiguous(get_ldz(), GRID_MPI_TYPE, &z_vector_);
  MPI_Type_commit(&z_vector_);

  MPI_Type_vector(get_ldy(), 1, get_ldz(), GRID_MPI_TYPE, &y_vector_);
  MPI_Type_commit(&y_vector_);
  
  MPI_Type_vector(get_ldx(), 1, get_ldy() * get_ldz(), 
                  GRID_MPI_TYPE, &x_vector_);
  MPI_Type_commit(&x_vector_);

  MPI_Type_contiguous(get_ldz() * get_ldy(), GRID_MPI_TYPE, &yz_plane_);
  MPI_Type_commit(&yz_plane_);

  MPI_Type_vector(get_ldx(), get_ldz(), get_ldy() * get_ldz(), 
                  GRID_MPI_TYPE, &xz_plane_);
  MPI_Type_commit(&xz_plane_);

  MPI_Type_hvector(get_ldx(), 1, sizeof(field_t) * get_ldz() * get_ldy(), 
                   y_vector_, &xy_plane_);
  MPI_Type_commit(&xy_plane_);

  types_alloced_ = true;
}


void Grid::alloc_grid()
{
  if (!define_)
  {
    cerr << "Unable to allocate grid data; the grid is not in define mode." << endl;
    return;
  }

  unsigned int sz = 0;

  sz = get_ldx() * get_ldy() * get_ldz();

  if (sz > 0) 
  {
    ex_ = new field_t[sz];
    ey_ = new field_t[sz];
    ez_ = new field_t[sz];
    
    hx_ = new field_t[sz];
    hy_ = new field_t[sz];
    hz_ = new field_t[sz];
    
    material_ = new unsigned int[sz];

    if (!ex_ || !ey_ || !ez_ || !hx_ || !hy_ || !hz_ || !material_) 
    {
      free_grid();
      throw MemoryException(); // Insufficient memory
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

  int num_mat = matlib.num_materials() + 1;
  Ca_ = new mat_coef_t[num_mat];
  Da_ = new mat_coef_t[num_mat];

  Cbx_ = new mat_coef_t[num_mat];
  Dbx_ = new mat_coef_t[num_mat];

  // Save some memory if possible. 
  if (get_deltay() == get_deltax())
  {
    Cby_ = Cbx_;
    Dby_ = Dbx_;
  } else {
    Cby_ = new mat_coef_t[num_mat];
    Dby_ = new mat_coef_t[num_mat];
  }

  if (get_deltaz() == get_deltax())
  {
    Cbz_ = Cbx_;
    Dbz_ = Dbx_;
  } else if (get_deltaz() == get_deltay()) {
    Cbz_ = Cby_;
    Dbz_ = Dby_;
  } else {
    Cbz_ = new mat_coef_t[num_mat];
    Dbz_ = new mat_coef_t[num_mat];
  }

  if (!Ca_ || !Da_ || !Cbx_ || !Cby_ || !Cbz_ || !Dbx_ || !Dby_ || !Dbz_) {
    free_material();
    throw MemoryException();
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
    mat_prop_t eps = (*iter).get_epsilon() * EPS_0;
    mat_prop_t sig = (*iter).get_sigma();
    mat_prop_t mu = (*iter).get_mu() * MU_0;
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
  
  if (!define_)
  {
    cerr << "Unable to define a box; the grid is not in define mode." << endl;
    return;
  }

  region_t r = global_to_local(x_start, x_stop,
                               y_start, y_stop,
                               z_start, z_stop);

  for (unsigned int i = r.xmin; i < r.xmax; i++)
  {
    for (unsigned int j = r.ymin; j < r.ymax; j++)
    {
      for (unsigned int k = r.zmin; k < r.zmax; k++)
      {
        material_[pi(i, j, k)] = mat_index;
      }
    }
  }
}

// Straight out of Taflove.
void Grid::update_e_field()
{
  if (define_)
  {
    cerr << "Unable to update fields; the grid is in define mode." << endl;
    return;
  }

  update_ex();
  update_ey();
  update_ez();
}

void Grid::update_h_field()
{
  if (define_)
  {
    cerr << "Unable to update fields; the grid is in define mode." << endl;
    return;
  }

  update_hx();
  update_hy();
  update_hz();

}

// Straight out of Taflove.
void Grid::update_ex() 
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

// Straight out of Taflove.
void Grid::update_ey() 
{
  unsigned int mid, idx;
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

        *ey = Ca_[mid] * *ey
          + Cbz_[mid] * (*hx - *(hx-1))
          + Cbx_[mid] * (*hz1 - *hz2);

        ey++;
        hx++;
        hz1++;
        hz2++;
        idx++;
      }
    }
  }

}

// Straight out of Taflove.
void Grid::update_ez() 
{
  unsigned int mid, idx;
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

        *ez = Ca_[mid] * *ez
          + Cbx_[mid] * (*hy1 - *hy2)
          + Cby_[mid] * (*hx1 - *hx2);

        ez++;
        hy1++; hy2++; hx1++; hx2++;
        idx++;
      }
    }
  }
}

// Straight out of Taflove.
void Grid::update_hx()
{
  unsigned int mid, idx;
  int i, j, k;
  field_t *hx, *ez1, *ez2, *ey;

  for (i = update_r_.xmin; i < update_r_.xmax; i++) {
    for (j = update_r_.ymin; j < update_r_.ymax - 1; j++) {

      idx = pi(i, j, update_r_.zmin);
      hx = &(hx_[idx]);
      ez1 = &(ez_[idx]);
      ez2 = &(ez_[pi(i, j+1, update_r_.zmin)]);
      ey = &(ey_[idx]);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (k = update_r_.zmin; k < update_r_.zmax - 1; k++) {
        mid = material_[idx];

        *hx = Da_[mid] * *hx
          + Dby_[mid] * (*ez1 - *ez2)
          + Dbz_[mid] * (*(ey+1) - *ey);

        hx++; idx++;
        ez1++; ez2++; ey++;
      }
    }
  }
}

// Straight out of Taflove.
void Grid::update_hy()
{
  unsigned int mid, idx;
  int i, j, k;
  field_t *hy, *ex, *ez1, *ez2;

  for (i = update_r_.xmin; i < update_r_.xmax - 1; i++) {
    for (j = update_r_.ymin; j < update_r_.ymax; j++) {

      idx = pi(i, j, update_r_.zmin);
      hy = &(hy_[idx]);
      ex = &(ex_[idx]);
      ez1 = &(ez_[pi(i+1, j, update_r_.zmin)]);
      ez2 = &(ez_[idx]);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (k = update_r_.zmin; k < update_r_.zmax - 1; k++) {
        mid = material_[idx];

        *hy = Da_[mid] * *hy
          + Dbz_[mid] * (*ex - *(ex + 1))
          + Dbx_[mid] * (*ez1 - *ez2);        
        
        hy++; idx++;
        ex++; ez1++; ez2++;
      }
    }
  }
}

// Straight out of Taflove.
void Grid::update_hz()
{
  unsigned int mid, idx;
  int i, j, k;
  field_t *hz1, *ey1, *ey2, *ex1, *ex2;

  for (i = update_r_.xmin; i < update_r_.xmax - 1; i++) {
    for (j = update_r_.ymin; j < update_r_.ymax - 1; j++) {

      idx = pi(i, j, update_r_.zmin);
      hz1 = &(hz_[idx]);
      ey1 = &(ey_[idx]);
      ey2 = &(ey_[pi(i+1, j, update_r_.zmin)]);
      ex1 = &(ex_[pi(i, j+1, update_r_.zmin)]);
      ex2 = &(ex_[idx]);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (k = update_r_.zmin; k < update_r_.zmax; k++) {
        mid = material_[idx];

        *hz1 = Da_[mid] * *hz1
          + Dbx_[mid] * (*ey1 - *ey2)
          + Dby_[mid] * (*ex1 - *ex2);

        hz1++; idx++;
        ey1++; ey2++;
        ex1++; ex2++;
      }
    }
  }
}

void Grid::apply_boundaries(FieldType type)
{
  if (define_)
  {
    cerr << "Unable to apply boundary conditions; the grid is in define mode." << endl;
    return;
  }

  info_.apply_boundaries(*this, type);
}


region_t Grid::global_to_local(region_t in) const
{
  region_t r;

  r.xmin = (get_lsx() > in.xmin) ? 0
    : in.xmin - get_lsx();
  r.ymin = (get_lsy() > in.ymin) ? 0
    : in.ymin - get_lsy();
  r.zmin = (get_lsz() > in.zmin) ? 0
    : in.zmin - get_lsz();

  r.xmax = (in.xmax >= get_lsx()) ? 
    ((in.xmax >= get_lsx() + get_ldx()) 
     ? get_ldx() : in.xmax - get_lsx() + 1)
    : 0;

  r.ymax = (in.ymax >= get_lsy()) ? 
    ((in.ymax >= get_lsy() + get_ldy()) 
     ? get_ldy() : in.ymax - get_lsy() + 1)
    : 0;

  r.zmax = (in.zmax >= get_lsz()) ? 
    ((in.zmax >= get_lsz() + get_ldz()) 
     ? get_ldz() : in.zmax - get_lsz() + 1)
    : 0;

  return r;
}

region_t Grid::global_to_local(unsigned int x_start, unsigned int x_stop, 
                             unsigned int y_start, unsigned int y_stop, 
                             unsigned int z_start, unsigned int z_stop) const
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

field_t *Grid::get_face_start(Face face, FieldComponent comp,
                              unsigned int offset)
{
  unsigned int idx = 0;
  field_t *ptr = 0;
  
  switch (face)
  {
  case FRONT: // x=dimx...
    idx = pi(get_ldx() - 1 - offset, 0, 0);
    break;
  case TOP: // z=dimz
    idx = pi(0, 0, get_ldz() - 1 - offset);
  case RIGHT: // y=dimy
    idx = pi(0, get_ldy() - 1 - offset, 0);
    break;
  case BACK: // x=0
    idx = pi(offset, 0, 0);
    break;
  case BOTTOM: // z=0
    idx = pi(0, 0, offset);
    break;
  case LEFT: // y=0
    idx = pi(0, offset, 0);
    break;
  }

  switch(comp)
  {
  case FC_EX:
    ptr = &(ex_[idx]);
    break;
  case FC_EY:
    ptr = &(ey_[idx]);
    break;
  case FC_EZ:
    ptr = &(ez_[idx]);
    break;
  case FC_HX:
    ptr = &(hx_[idx]);
    break;
  case FC_HY:
    ptr = &(hy_[idx]);
    break;
  case FC_HZ:
    ptr = &(hz_[idx]);
    break;
  }

  return ptr;
}

field_t *Grid::get_face_start(Face face, FieldComponent comp,
                              point_t p) 
{
  unsigned int idx = 0;
  field_t *ptr = 0;

  switch (face)
  {
  case FRONT:
  case BACK:
    idx = pi(p.x, 0, 0);    
    break;

  case TOP:
  case BOTTOM:
    idx = pi(0, 0, p.z);    
    break;

  case LEFT:
  case RIGHT:
    idx = pi(0, p.y, 0);
    break;
  }

  switch(comp)
  {
  case FC_EX:
    ptr = &(ex_[idx]);
    break;
  case FC_EY:
    ptr = &(ey_[idx]);
    break;
  case FC_EZ:
    ptr = &(ez_[idx]);
    break;
  case FC_HX:
    ptr = &(hx_[idx]);
    break;
  case FC_HY:
    ptr = &(hy_[idx]);
    break;
  case FC_HZ:
    ptr = &(hz_[idx]);
    break;
  }

  return ptr;
}

MPI_Datatype Grid::get_plane_dt(Face face) const
{
  MPI_Datatype t;

  switch (face)
  {
  case FRONT:
  case BACK:
    t = yz_plane_;
    break;

  case TOP:
  case BOTTOM:
    t = xy_plane_;
    break;

  case LEFT:
  case RIGHT:
    t = xz_plane_;
    break;
  }

  return t;
}

const field_t *Grid::get_pointer(point_t point, 
                                 FieldComponent field_comp) const
{
  field_t *ret = 0;
  unsigned int idx = pi(point.x, point.y, point.z);

  switch(field_comp)
  {
  case FC_EX:
    if (ex_)
      ret = &(ex_[idx]);
    break;
  case FC_EY:
    if (ey_)
      ret = &(ey_[idx]);
    break;
  case FC_EZ:
    if (ez_)
      ret = &(ez_[idx]);
    break;
  case FC_HX:
    if (hx_)
      ret = &(hx_[idx]);
    break;
  case FC_HY:
    if (hy_)
      ret = &(hy_[idx]);
    break;
  case FC_HZ:
    if (hz_)
      ret = &(hz_[idx]);
    break;
  }

  return ret;
}
