/* 
   Phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004-2005 Matt Hughes <mhughe@uvic.ca>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#include "config.h"

#include "Grid.hh"
#include "Globals.hh"

#include "Boundaries/Ewall.hh"
#include "Boundaries/Hwall.hh"
#include "Constants.hh"
#include "Exceptions.hh"

#include <mpi.h>
#include <cmath>
#include <cstring> // for memset

#ifdef USE_OPENMP
#include <omp.h>
#endif

Grid::Grid() 
  : num_materials_(0),
#ifdef OLD_MATERIAL_DATA
    Ca_(0), Cbx_(0), Cby_(0), Cbz_(0),
    Da_(0), Dbx_(0), Dby_(0), Dbz_(0),
#else
    C_(0), D_(0),
#endif
    ex_(0), ey_(0), ez_(0), hx_(0), hy_(0), hz_(0), 
    material_(0), types_alloced_(false), define_(true)
    //geometries_(0), num_geoms_(0)
{

}

Grid::Grid(const Grid &rhs)
{
  *this = rhs;
}

const Grid &Grid::operator=(const Grid &rhs)
{
  info_ = rhs.info_;

  update_ex_r_ = rhs.update_ex_r_;
  update_ey_r_ = rhs.update_ey_r_;
  update_ez_r_ = rhs.update_ez_r_;

  update_hx_r_ = rhs.update_hx_r_;
  update_hy_r_ = rhs.update_hy_r_;
  update_hz_r_ = rhs.update_hz_r_;

  num_materials_ = rhs.num_materials_;

#ifdef OLD_MATERIAL_DATA
  Ca_ = Cbx_ = Cby_ = Cbz_ = Da_ = Dbx_ = Dby_ = Dbz_ = 0;
#else
  C_ = D_ = 0;
#endif

  ex_ = ey_ = ez_ = hx_ = hy_ = hz_ = 0;
  material_lib_ = rhs.material_lib_;
  material_ = 0;
  xy_plane_ = rhs.xy_plane_;
  yz_plane_ = rhs.yz_plane_;
  xz_plane_ = rhs.xz_plane_;
  x_vector_ = rhs.x_vector_;
  y_vector_ = rhs.y_vector_;
  z_vector_ = rhs.z_vector_;
  define_ = rhs.define_;

//   geometries_ = 0;
//   num_geoms_ = 0;

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
    float temp = C * sqrt( 1/(get_deltax() * get_deltax()) + 
                           1/(get_deltay() * get_deltay()) + 
                           1/(get_deltaz() * get_deltaz()));

    if (get_deltat() > 1/temp)
      throw StabilityException();
    
    if (!setup_only)
      alloc_grid();

    // Set up the grid geometry. 
    if (!setup_only && pg_)
    {
      info_.default_mat_ = pg_->get_grid_material_id();

      point size = (*pg_).get_grid_size();
      point centre = (*pg_).get_grid_centre();

      float x = centre.x - size.x / 2 + get_deltax() / 2 
        + get_deltax() * info_.start_x_;
      float y = centre.y - size.y / 2 + get_deltay() / 2 
        + get_deltay() * info_.start_y_;
      float z = centre.z - size.z / 2 + get_deltaz() / 2 
        + get_deltaz() * info_.start_z_;

      for (unsigned int i = 0; i < info_.dimx_; i++)
      {
        y = centre.y - size.y / 2 + get_deltay() * info_.start_y_
          + get_deltay() / 2 ;
        for (unsigned int j = 0; j < info_.dimy_; j++)
        {
          z = centre.z - size.z / 2 + get_deltaz() * info_.start_z_
            + get_deltaz() / 2 ;
          for (unsigned int k = 0; k < info_.dimz_; k++)
          {
            material_[pi(i,j,k)] = pg_->get_material_id(x, y, z);

            z += info_.deltax_;
          }
          y += info_.deltay_;
        }
        x += info_.deltax_;
      }

      // SET UP Data structures for dealing with blocks of memory for different
      // kinds of dispersion relationships. 

      //   num_geoms_ = geoms.size();
      //   vector<Geometry *>::iterator iter;
      //   vector<Geometry *>::iterator iter_e = geoms.end();
      //   int idx = 0;
  
      //   if (geometries_)
      //     delete geometries_;

      //   geometries_ = new Geometry*[num_geoms_];

      //   for (iter = geoms.begin(), idx = 0; iter != iter_e; ++iter, idx++)
      //     geometries_[idx] = *iter;
    }

    init_datatypes();

    // Calculate update region_t by considering the thickness of the PML's. 
    region_t update_r;
    update_r.xmin = 0;
    update_r.xmax = info_.dimx_;
    update_r.ymin = 0;
    update_r.ymax = info_.dimy_;
    update_r.zmin = 0;
    update_r.zmax = info_.dimz_;

    update_ex_r_ = update_ey_r_ = update_ez_r_ = update_r;
    update_hx_r_ = update_hy_r_ = update_hz_r_ = update_r;

    update_ex_r_.ymin++;
    update_ex_r_.zmin++;

    update_ey_r_.xmin++;
    update_ey_r_.zmin++;

    update_ez_r_.xmin++;
    update_ez_r_.ymin++;

    update_hx_r_.ymax--;
    update_hx_r_.zmax--;

    update_hy_r_.xmax--;
    update_hy_r_.zmax--;

    update_hz_r_.xmax--;
    update_hz_r_.ymax--;

    // The domain decomp algorithm will take care of assigning the
    // boundary conditions sensibly, so we don't have to worry about
    // wether or not we are really on a boundary that has thickness
    // (i.e. a PML face)
    unsigned int thickness = 0;
    for (int i = 0; i < 6; i++)
    {
      thickness = info_.get_face_thickness(static_cast<Face>(i));

      if (thickness > 0) 
      {
        switch (i) {
        case FRONT:
          update_ex_r_.xmax -= (thickness + 1);
          update_ey_r_.xmax -= (thickness + 1);
          update_ez_r_.xmax -= (thickness + 1);

          update_hx_r_.xmax -= (thickness + 1);
          update_hy_r_.xmax -= thickness;
          update_hz_r_.xmax -= thickness;
          break;

        case BACK:
          update_ex_r_.xmin += thickness;
          update_ey_r_.xmin += thickness;
          update_ez_r_.xmin += thickness;

          update_hx_r_.xmin += (thickness + 1);
          update_hy_r_.xmin += thickness;
          update_hz_r_.xmin += thickness;
          break;

        case TOP:
          update_ex_r_.zmax -= (thickness + 1);
          update_ey_r_.zmax -= (thickness + 1);
          update_ez_r_.zmax -= (thickness + 1);

          update_hx_r_.zmax -= thickness;
          update_hy_r_.zmax -= thickness;
          update_hz_r_.zmax -= (thickness + 1);
          break;

        case BOTTOM:
          update_ex_r_.zmin += thickness;
          update_ey_r_.zmin += thickness;
          update_ez_r_.zmin += thickness;

          update_hx_r_.zmin += thickness;
          update_hy_r_.zmin += thickness;
          update_hz_r_.zmin += (thickness + 1);
          break;

        case LEFT:
          update_ex_r_.ymin += thickness;
          update_ey_r_.ymin += thickness;
          update_ez_r_.ymin += thickness;

          update_hx_r_.ymin += thickness;
          update_hy_r_.ymin += (thickness + 1);
          update_hz_r_.ymin += thickness;
          break;

        case RIGHT:
          update_ex_r_.ymax -= (thickness + 1);
          update_ey_r_.ymax -= (thickness + 1);
          update_ez_r_.ymax -= (thickness + 1);

          update_hx_r_.ymax -= thickness;
          update_hy_r_.ymax -= (thickness + 1);
          update_hz_r_.ymax -= thickness;
          break;
        }
      }
      // Initialize the Boundary conditions. 
      BoundaryCond *p = &info_.get_boundary(static_cast<Face>(i));

      if (p) {
        p->init(*this, static_cast<Face>(i));

        // Check ajacent faces and see if there are any subdomains
        // that need to have data shared across them.

//         cout << "RANK " << MPI_RANK << ", "  
//              << face_string(static_cast<Face>(i)) 
//              << " shares data across [";

        for (int j = 0; j < 6; j++)
        {
          if (j != i 
              && ((j % 2 && (j-1) != i) // j is odd
                  || ( !(j % 2) && (j+1) != i))) // j is even
          {
            if (info_.get_bc_type(static_cast<Face>(j)) == SUBDOMAIN)
            {
              SubdomainBc *sd = dynamic_cast<SubdomainBc *>
                (&info_.get_boundary(static_cast<Face>(j)));

//               cout << face_string(static_cast<Face>(j)) << ", ";

              p->add_sd_bcs(sd, static_cast<Face>(i), 
                            static_cast<Face>(j));
            }
          }
        }

//         cout << "\b\b]" << endl;
      }

      // Tell subdomains about Grid data that needs to be exchanged
      SubdomainBc *sd = 
        dynamic_cast<SubdomainBc *>(&info_.get_boundary(static_cast<Face>(i)));
      if (sd)
      {
        setup_subdomain_data(sd, static_cast<Face>(i));
      }
    }

    // Init subdomains
    for (int i = 0; i < 6; i++)
    {
      SubdomainBc *sd = 
        dynamic_cast<SubdomainBc *>(&info_.get_boundary(static_cast<Face>(i)));
      
      if(sd)
      {
        sd->sd_init(*this, static_cast<Face>(i));
      }
    }

#ifdef DEBUG    
  cout << "Grid Update region:"
       << "\n\tEx, x: " << update_ex_r_.xmin << " -> " 
       << update_ex_r_.xmax
       << ", y: " << update_ex_r_.ymin << " -> " 
       << update_ex_r_.ymax
       << ", z: " << update_ex_r_.zmin << " -> " 
       << update_ex_r_.zmax
       << "\n\tEy, x: " << update_ey_r_.xmin << " -> " 
       << update_ey_r_.xmax
       << ", y: " << update_ey_r_.ymin << " -> " 
       << update_ey_r_.ymax
       << ", z: " << update_ey_r_.zmin << " -> " 
       << update_ey_r_.zmax
       << "\n\tEz, x: " << update_ez_r_.xmin << " -> " 
       << update_ez_r_.xmax
       << ", y: " << update_ez_r_.ymin << " -> " 
       << update_ez_r_.ymax
       << ", z: " << update_ez_r_.zmin << " -> " 
       << update_ez_r_.zmax 
       << "\n\tHx, x: " << update_hx_r_.xmin << " -> " 
       << update_hx_r_.xmax
       << ", y: " << update_hx_r_.ymin << " -> " 
       << update_hx_r_.ymax
       << ", z: " << update_hx_r_.zmin << " -> " 
       << update_hx_r_.zmax
       << "\n\tHy, x: " << update_hy_r_.xmin << " -> " 
       << update_hy_r_.xmax
       << ", y: " << update_hy_r_.ymin << " -> " 
       << update_hy_r_.ymax
       << ", z: " << update_hy_r_.zmin << " -> " 
       << update_hy_r_.zmax
       << "\n\tHz, x: " << update_hz_r_.xmin << " -> " 
       << update_hz_r_.xmax
       << ", y: " << update_hz_r_.ymin << " -> " 
       << update_hz_r_.ymax
       << ", z: " << update_hz_r_.zmin << " -> " 
       << update_hz_r_.zmax << endl;
#endif

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

#ifdef OLD_MATERIAL_DATA
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
#else
  if (C_)
    delete[] C_;

  if (D_)
    delete[] D_;
  
  C_ = D_ = 0;
#endif
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

  MPI_Type_contiguous(info_.dimz_, GRID_MPI_TYPE, &z_vector_);
  MPI_Type_commit(&z_vector_);

  MPI_Type_vector(info_.dimy_, 1, info_.dimz_, GRID_MPI_TYPE, &y_vector_);
  MPI_Type_commit(&y_vector_);
  
  MPI_Type_vector(info_.dimx_, 1, info_.dimy_ * info_.dimz_, 
                  GRID_MPI_TYPE, &x_vector_);
  MPI_Type_commit(&x_vector_);

  MPI_Type_contiguous(info_.dimz_ * info_.dimy_, GRID_MPI_TYPE, &yz_plane_);
  MPI_Type_commit(&yz_plane_);

  MPI_Type_vector(info_.dimx_, info_.dimz_, info_.dimy_ * info_.dimz_, 
                  GRID_MPI_TYPE, &xz_plane_);
  MPI_Type_commit(&xz_plane_);

  MPI_Type_hvector(info_.dimx_, 1, sizeof(field_t) * info_.dimz_ * info_.dimy_, 
                   y_vector_, &xy_plane_);
  MPI_Type_commit(&xy_plane_);

#ifdef DEBUG
  int sz; 
  MPI_Aint extnt;

  MPI_Type_size(x_vector_, &sz);
  MPI_Type_extent(x_vector_, &extnt);

  cout << "Grid::init_datatypes(): x vector is " << sz << " bytes in size, "
       << extnt << " in extent. " << endl;

  MPI_Type_size(y_vector_, &sz);
  MPI_Type_extent(y_vector_, &extnt);

  cout << "Grid::init_datatypes(): y vector is " << sz << " bytes in size, "
       << extnt << " in extent. " << endl;

  MPI_Type_size(z_vector_, &sz);
  MPI_Type_extent(z_vector_, &extnt);

  cout << "Grid::init_datatypes(): z vector is " << sz << " bytes in size, "
       << extnt << " in extent. " << endl;

  MPI_Type_size(yz_plane_, &sz);
  MPI_Type_extent(yz_plane_, &extnt);

  cout << "Grid::init_datatypes(): yz plane is " << sz << " bytes in size, "
       << extnt << " in extent. " << endl;

  MPI_Type_size(xz_plane_, &sz);
  MPI_Type_extent(xz_plane_, &extnt);

  cout << "Grid::init_datatypes(): xz plane is " << sz << " bytes in size, "
       << extnt << " in extent. " << endl;

  MPI_Type_size(xy_plane_, &sz);
  MPI_Type_extent(xy_plane_, &extnt);

  cout << "Grid::init_datatypes(): xy plane is " << sz << " bytes in size, "
       << extnt << " in extent. " << endl;
#endif

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

  sz = info_.dimx_ * info_.dimy_ * info_.dimz_;

  if (sz > 0) 
  {
    ex_ = new field_t[sz];
    ey_ = new field_t[sz];
    ez_ = new field_t[sz];
    
    hx_ = new field_t[sz];
    hy_ = new field_t[sz];
    hz_ = new field_t[sz];
    
    material_ = new mat_idx_t[sz];

    if (!ex_ || !ey_ || !ez_ || !hx_ || !hy_ || !hz_ || !material_) 
    {
      free_grid();
      throw MemoryException(); // Insufficient memory
    }
    
    memset(ex_, 0, sizeof(field_t) * sz);
    memset(ey_, 0, sizeof(field_t) * sz);
    memset(ez_, 0, sizeof(field_t) * sz);

    memset(hx_, 0, sizeof(field_t) * sz);
    memset(hy_, 0, sizeof(field_t) * sz);
    memset(hz_, 0, sizeof(field_t) * sz);

    memset(material_, 0, sizeof(mat_idx_t) * sz);

#ifdef DEBUG
    cout << "Sucessfully allocated " 
	 << sz * sizeof(field_t) * 6 + sz * sizeof(unsigned int) 
      + sz * sizeof(mat_idx_t)
	 << " bytes of memory for use by the grid." << endl;
#endif
  }
}

void Grid::load_geometry(const ProblemGeometry *pg)
{
  // Loop through the voxels and ask the problem geometry what
  // material id should be used at each. THIS MUST BE CHANGED TO
  // SUPPORT GRADED MESHES.

  pg_ = pg;

}

void Grid::load_materials(shared_ptr<MaterialLib> matlib)
{
  if (!define_)
  {
    cerr << "Unable to load material data; the grid is not in define mode." << endl;
    return;
  }
  
  material_lib_ = matlib;
  
  // Clear up any material data that may already be loaded
  free_material();

  int num_mat = (*matlib).num_materials() + 1;

#ifdef OLD_MATERIAL_DATA
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

  memset(Ca_, 0, sizeof(mat_coef_t) * num_mat);
  memset(Da_, 0, sizeof(mat_coef_t) * num_mat);
  memset(Cbx_, 0, sizeof(mat_coef_t) * num_mat);
  memset(Cbz_, 0, sizeof(mat_coef_t) * num_mat);
  memset(Cby_, 0, sizeof(mat_coef_t) * num_mat);
  memset(Dbx_, 0, sizeof(mat_coef_t) * num_mat);
  memset(Dby_, 0, sizeof(mat_coef_t) * num_mat);
  memset(Dbz_, 0, sizeof(mat_coef_t) * num_mat);
#else
  C_ = new mat_coef_t[num_mat * 4];
  D_ = new mat_coef_t[num_mat * 4];

  memset(C_, 0, sizeof(mat_coef_t) * num_mat * 4);
  memset(D_, 0, sizeof(mat_coef_t) * num_mat * 4);
#endif
  map<string, Material>::iterator iter = (*matlib).materials_.begin();
  map<string, Material>::iterator iter_e = (*matlib).materials_.end();

  // The first one is always PEC
  unsigned int index = 0;

//   Ca_[index] = 1;
//   Cbx_[index] = Cby_[index] = Cbz_[index] = 0;

//   Da_[index] = 1;
//   Dbx_[index] = Dby_[index] = Dbz_[index] = 0;
  
//   ++index;
  while (iter != iter_e) 
  {
    // Set the index so that the geometry objects know about it. 
    (*iter).second.set_id(index);
    
    mat_prop_t eps = ((*iter).second).get_epsilon() * EPS_0;
    mat_prop_t sig = ((*iter).second).get_sigma();
    mat_prop_t mu = ((*iter).second).get_mu() * MU_0;
    mat_prop_t sigs = ((*iter).second).get_sigma_star();

    if (((*iter).second).is_pec())
    {
#ifdef OLD_MATERIAL_DATA
      Ca_[index] = 1;
      Cbx_[index] = Cby_[index] = Cbz_[index] = 0;

      Da_[index] = 1;
      Dbx_[index] = Dby_[index] = Dbz_[index] = 0;
#else
      C_[index * 4] = 1;
      C_[index * 4 + 1] = 0;
      C_[index * 4 + 2] = 0;
      C_[index * 4 + 3] = 0;

      D_[index * 4] = 1;
      D_[index * 4 + 1] = 0;
      D_[index * 4 + 2] = 0;
      D_[index * 4 + 3] = 0;
#endif
    } 
    else if (eps == 0 || mu == 0)
    {
      cerr << "Something is wrong with the material library:\n" 
           << " -> Material cannot have permittivities or permeabilities\n"
           << "    of zero. Perfect electric conductor can have eps=0, \n"
           << "    but that is a special material defined by Phred.\n\n"
           << "Program aborting. Check material library." << endl;
      MPI_Abort(MPI_COMM_PHRED, 1);
    }
    else 
    {
#ifdef OLD_MATERIAL_DATA
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
#else
      C_[index * 4] = (1 - (sig * get_deltat() * 0.5)/eps) / 
                      (1 + (sig * get_deltat() * 0.5)/eps);

      D_[index * 4] = (1 - (sigs * get_deltat() * 0.5)/mu) / 
                      (1 + (sigs * get_deltat() * 0.5)/mu);

      C_[index * 4 + 1] = (get_deltat() / (eps * get_deltax())) / 
                          (1 + (sig * get_deltat() * 0.5)/eps);

      D_[index * 4 + 1] = (get_deltat() / (mu * get_deltax())) / 
                          (1 + (sigs * get_deltat() * 0.5)/mu);

      C_[index * 4 + 2] = (get_deltat() / (eps * get_deltay())) / 
                          (1 + (sig * get_deltat() * 0.5)/eps);

      D_[index * 4 + 2] = (get_deltat() / (mu * get_deltay())) / 
                          (1 + (sigs * get_deltat() * 0.5)/mu);

      C_[index * 4 + 3] = (get_deltat() / (eps * get_deltaz())) / 
                          (1 + (sig * get_deltat() * 0.5)/eps);

      D_[index * 4 + 3] = (get_deltat() / (mu * get_deltaz())) / 
                          (1 + (sigs * get_deltat() * 0.5)/mu);      
#endif
    }

#ifdef DEBUG
//     cerr << "Material '" << (*iter).second.get_name() << "', index: " 
//          << index << "\n\tCa_ = " << Ca_[index]
//          << ", Da_ = " << Da_[index] 
//          << "\n\tCbx_ = " << Cbx_[index] << ", Cby_" << Cby_[index]
//          << ", Cbz_ = " << Cbz_[index]
//          << "\n\tDbx_ = " << Dbx_[index] << ", Dby_" << Dby_[index]
//          << ", Dbz_ = " << Dbz_[index]
//          << endl;
#endif

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


// Straight out of Taflove.
void Grid::update_e_field()
{
  if (setup_only)
    return;

  if (define_)
  {
    cerr << "Unable to update fields; the grid is in define mode." << endl;
    return;
  }

  update_ex(update_ex_r_);
  update_ey(update_ey_r_);
  update_ez(update_ez_r_);
}

void Grid::update_h_field()
{
  if (setup_only)
    return;

  if (define_)
  {
    cerr << "Unable to update fields; the grid is in define mode." << endl;
    return;
  }

  update_hx(update_hx_r_);
  update_hy(update_hy_r_);
  update_hz(update_hz_r_);

}

// Straight out of Taflove.
void Grid::update_ex(region_t update_r) 
{
  unsigned int idx, idx2;
  mat_idx_t mid;
  field_t *ex, *hz1, *hz2, *hy;

  // Inner part
#ifdef USE_OPENMP
#pragma omp parallel  private(mid, idx, idx2, ex, hz1, hz2, hy)
#endif
  {
#ifdef USE_OPENMP
#pragma omp for
#endif
    for (loop_idx_t i = update_r.xmin; i < update_r.xmax; i++) {
      for (loop_idx_t j = update_r.ymin; j < update_r.ymax; j++) {
        
        idx = pi(i, j, update_r.zmin);
        idx2 = pi(i, j-1, update_r.zmin);

        ex = &(ex_[idx]);
        hz1 = &(hz_[idx]);
        hz2 = &(hz_[idx2]);
        hy = &(hy_[idx]);
        
        for (loop_idx_t k = update_r.zmin; k < update_r.zmax; k++) {
          mid = material_[idx];
          
          *ex = get_Ca(mid) * *ex
            + get_Cby(mid) * (*hz1 - *hz2)
            + get_Cbz(mid) * (*(hy - 1) - *hy);
          
          ex++;
          hz1++;
          hz2++;
          hy++;
          idx++;
        }
      }
    }
  }

}

// Straight out of Taflove.
void Grid::update_ey(region_t update_r) 
{
  unsigned int idx;
  mat_idx_t mid;
  field_t *ey, *hx, *hz1, *hz2;

  // Inner part
#ifdef USE_OPENMP
#pragma omp parallel private(mid, idx, ey, hx, hz1, hz2)
#endif
  {
#ifdef USE_OPENMP
#pragma omp for
#endif
    for (loop_idx_t i = update_r.xmin; i < update_r.xmax; i++) {
      for (loop_idx_t j = update_r.ymin; j < update_r.ymax; j++) {
        
        idx = pi(i, j, update_r.zmin);
        hz1 = &(hz_[pi(i-1, j, update_r.zmin)]);

        hz2 = &(hz_[idx]);
        ey = &(ey_[idx]);
        hx = &(hx_[idx]);

        for (loop_idx_t k = update_r.zmin; k < update_r.zmax; k++) {
          mid = material_[idx];
          
          *ey = get_Ca(mid) * *ey
            + get_Cbz(mid) * (*hx - *(hx-1))
            + get_Cbx(mid) * (*hz1 - *hz2);
          
          ey++;
          hx++;
          hz1++;
          hz2++;
          idx++;
        }
      }
    }
  }
}

// Straight out of Taflove.
void Grid::update_ez(region_t update_r) 
{
  unsigned int idx;
  mat_idx_t mid;
  field_t *ez, *hy1, *hy2, *hx1, *hx2;
  
  // Inner part
#ifdef USE_OPENMP
#pragma omp parallel private(mid, idx, ez, hy1, hy2, hx1, hx2)
#endif
  {
#ifdef USE_OPENMP
#pragma omp for
#endif
    for (loop_idx_t i = update_r.xmin; i < update_r.xmax; i++) {
      for (loop_idx_t j = update_r.ymin; j < update_r.ymax; j++) {

        idx = pi(i, j, update_r.zmin);
        hy2 = &(hy_[pi(i-1, j, update_r.zmin)]);
        hx1 = &(hx_[pi(i, j-1, update_r.zmin)]);

        hx2 = &(hx_[idx]);
        ez = &(ez_[idx]);
        hy1 = &(hy_[idx]);

        for (loop_idx_t k = update_r.zmin; k < update_r.zmax; k++) {
          mid = material_[idx];
          
          *ez = get_Ca(mid) * *ez
            + get_Cbx(mid) * (*hy1 - *hy2)
            + get_Cby(mid) * (*hx1 - *hx2);

          ez++;
          hy1++; hy2++; hx1++; hx2++;
          idx++;
        }
      }
    }
  }
}

// Straight out of Taflove.
void Grid::update_hx(region_t update_r)
{
  unsigned int idx;
  mat_idx_t mid;
  field_t *hx, *ez1, *ez2, *ey;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, idx, hx, ez1, ez2, ey)
#endif
  {
#ifdef USE_OPENMP
#pragma omp for
#endif
    for (loop_idx_t i = update_r.xmin; i < update_r.xmax; i++) {
      for (loop_idx_t j = update_r.ymin; j < update_r.ymax; j++) {
        
        idx = pi(i, j, update_r.zmin);
        ez2 = &(ez_[pi(i, j+1, update_r.zmin)]);

        ey = &(ey_[idx]);
        hx = &(hx_[idx]);
        ez1 = &(ez_[idx]);

        for (loop_idx_t k = update_r.zmin; k < update_r.zmax; k++) {
          mid = material_[idx];
          
          *hx = get_Da(mid) * *hx
            + get_Dby(mid) * (*ez1 - *ez2)
            + get_Dbz(mid) * (*(ey+1) - *ey);
          
          hx++; idx++;
          ez1++; ez2++; ey++;
        }
      }
    }
  }

}

// Straight out of Taflove.
void Grid::update_hy(region_t update_r)
{
  unsigned int idx;
  mat_idx_t mid;
  field_t *hy, *ex, *ez1, *ez2;


#ifdef USE_OPENMP
#pragma omp parallel private(mid, idx, hy, ex, ez1, ez2)
#endif
  {
#ifdef USE_OPENMP
#pragma omp for
#endif
    for (loop_idx_t i = update_r.xmin; i < update_r.xmax; i++) {
      for (loop_idx_t j = update_r.ymin; j < update_r.ymax; j++) {

        idx = pi(i, j, update_r.zmin);
        ez1 = &(ez_[pi(i+1, j, update_r.zmin)]);

        hy = &(hy_[idx]);
        ex = &(ex_[idx]);
        ez2 = &(ez_[idx]);

        for (loop_idx_t k = update_r.zmin; k < update_r.zmax; k++) {
          mid = material_[idx];
          
          *hy = get_Da(mid) * *hy
            + get_Dbz(mid) * (*ex - *(ex + 1))
            + get_Dbx(mid) * (*ez1 - *ez2);        
          
          hy++; idx++;
          ex++; ez1++; ez2++;
        }
      }
    }
  }
}

// Straight out of Taflove.
void Grid::update_hz(region_t update_r)
{
  unsigned int idx;
  mat_idx_t mid;

  field_t *hz1, *ey1, *ey2, *ex1, *ex2;

#ifdef USE_OPENMP
#pragma omp parallel private(mid, idx, hz1, ey1, ey2, ex1, ex2)
#endif
  {
#ifdef USE_OPENMP
#pragma omp for
#endif
    for (loop_idx_t i = update_r.xmin; i < update_r.xmax; i++) {
      for (loop_idx_t j = update_r.ymin; j < update_r.ymax; j++) {
        
        idx = pi(i, j, update_r.zmin);
        ey2 = &(ey_[pi(i+1, j, update_r.zmin)]);
        ex1 = &(ex_[pi(i, j+1, update_r.zmin)]);

        ex2 = &(ex_[idx]);
        hz1 = &(hz_[idx]);
        ey1 = &(ey_[idx]);

        for (loop_idx_t k = update_r.zmin; k < update_r.zmax; k++) {
          mid = material_[idx];
          
          *hz1 = get_Da(mid) * *hz1
            + get_Dbx(mid) * (*ey1 - *ey2)
            + get_Dby(mid) * (*ex1 - *ex2);
          
          hz1++; idx++;
          ey1++; ey2++;
          ex1++; ex2++;
        }
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

  if (!setup_only)
    info_.apply_boundaries(*this, type);
}

field_t *Grid::get_face_start(Face face, FieldComponent comp,
                              unsigned int offset) const
{
  unsigned int idx = 0;
  field_t *ptr = 0;
  
  switch (face)
  {
  case FRONT: // x=dimx...
    idx = pi(info_.dimx_ - 1 - offset, 0, 0);
    break;
  case TOP: // z=dimz
    idx = pi(0, 0, info_.dimz_ - 1 - offset);
    break;
  case RIGHT: // y=dimy
    idx = pi(0, info_.dimy_ - 1 - offset, 0);
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

const mat_idx_t *Grid::get_material_ptr(grid_point point) const
{
  return &(material_[pi(point.x, point.y, point.z)]);
}

const field_t *Grid::get_pointer(grid_point point, 
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

void Grid::setup_subdomain_data(SubdomainBc *sd, Face face)
{
  RxTxData rxtx;
  MPI_Datatype t = get_plane_dt(face);
  rxtx.set_field_type(E);
  
  rxtx.set_datatype(t);
  rxtx.set_tx_ptr(get_face_start(face, FC_EX, 1));
  rxtx.set_rx_ptr(get_face_start(face, FC_EX, 0));
  sd->add_tx_rx_data(rxtx);
  
  rxtx.set_tx_ptr(get_face_start(face, FC_EY, 1));
  rxtx.set_rx_ptr(get_face_start(face, FC_EY, 0));
  sd->add_tx_rx_data(rxtx);
  
  rxtx.set_tx_ptr(get_face_start(face, FC_EZ, 1));
  rxtx.set_rx_ptr(get_face_start(face, FC_EZ, 0));
  sd->add_tx_rx_data(rxtx);
  
  rxtx.set_field_type(H);
  
  rxtx.set_tx_ptr(get_face_start(face, FC_HX, 1));
  rxtx.set_rx_ptr(get_face_start(face, FC_HX, 0));
  sd->add_tx_rx_data(rxtx);
  
  rxtx.set_tx_ptr(get_face_start(face, FC_HY, 1));
  rxtx.set_rx_ptr(get_face_start(face, FC_HY, 0));
  sd->add_tx_rx_data(rxtx);

  rxtx.set_tx_ptr(get_face_start(face, FC_HZ, 1));
  rxtx.set_rx_ptr(get_face_start(face, FC_HZ, 0));
  sd->add_tx_rx_data(rxtx);

#ifdef DEBUG
  cout << "Grid::setup_subdomain_data(): subdomain bc on " 
       << face_string(face) << " on rank " << sd->get_rank() 
       << " talking with " << sd->get_neighbour() << endl;
#endif
}

grid_point Grid::get_global_cell(point p) const
{
  return get_global_cell(p.x, p.y, p.z);
}

// MCH, 2005-02-08: Modified to use global lengths rather than local
// when checking if a point is outside of the grid.
grid_point Grid::get_global_cell(float x, float y, float z) const
{
  grid_point ret; 

  if (pg_)
  {
    point size = (*pg_).get_grid_size();
    point centre = (*pg_).get_grid_centre();

    // Lower left back corner of the grid
    float xs = centre.x - size.x / 2;
    float ys = centre.y - size.y / 2;
    float zs = centre.z - size.z / 2;
    
    unsigned int i = 0; 
    unsigned int j = 0; 
    unsigned int k = 0;

    // Will need modification for graded meshes. 
    if (x > xs)
    {
      i = static_cast<unsigned int>(floor((x - xs) / get_deltax()));
      if (i >= get_gdx())
        i = get_gdx() - 1;
      //if (i > get_ldx_sd())
      //  i = get_ldx_sd() - 1;
    }

    if (y > ys)
    {
      j = static_cast<unsigned int>(floor((y - ys) / get_deltay()));
      if (j >= get_gdy())
        j = get_gdy() - 1;
      //if (j > get_ldy_sd())
      //  j = get_ldy_sd() - 1;
    }

    if (z > zs)
    {
      k = static_cast<unsigned int>(floor((z - zs) / get_deltaz()));
      if (k >= get_gdz())
        k = get_gdz() - 1;
      //if (k > get_ldz_sd())
      //  k = get_ldz_sd() - 1;
    }

    ret.x = i;
    ret.y = j; 
    ret.z = k;

  } else {
    cerr << "Grid::get_global_cell: NO PROBLEM GEOMETRY IS LOADED INTO THE GRID!" << endl;
  }

  return ret;
}

shared_ptr<CellSet> Grid::get_cellset(const CSGBox &box) const
{
  shared_ptr<CellSet> cells = shared_ptr<CellSet>(new CellSet());
  
  // Fill out the global set
  point centre = box.get_centre();
  point size = box.get_size();
  shared_ptr<Block> global = cells->global_;

  float xs = centre.x - size.x / 2;
  float ys = centre.y - size.y / 2;
  float zs = centre.z - size.z / 2;
  
  float xe = centre.x + size.x / 2;
  float ye = centre.y + size.y / 2;
  float ze = centre.z + size.z / 2;
  
  grid_point start = get_global_cell(xs, ys, zs);
  grid_point end = get_global_cell(xe, ye, ze);
  
  global->xmin_ = start.x; global->xmax_ = end.x;
  global->ymin_ = start.y; global->ymax_ = end.y;
  global->zmin_ = start.z; global->zmax_ = end.z;

  global->xlen_ = global->xmax_ - global->xmin_ + 1;
  global->ylen_ = global->ymax_ - global->ymin_ + 1;
  global->zlen_ = global->zmax_ - global->zmin_ + 1;

  if (global->xlen_ == 0 || global->ylen_ == 0 || global->zlen_ == 0)
  {
    global->has_data_ = false;

    for (int i = 0; i < 6; i++)
      global->faces_[i] = false;
  }

  // Local, with ghosts
  cells->local_ghost_ = global_to_local_ghost(global);

  // Local, without ghosts
  cells->local_ = global_to_local(global);

// #ifdef DEBUG
//   cerr << "Grid::get_cellset, global block: \n" 
//        << *global << "\nlocal block, no ghosts:\n"
//        << *(cells->local_) << "\nlocal block, with ghosts:\n"
//        << *(cells->local_ghost_) << endl << endl;
// #endif

  return cells;
}

shared_ptr<Block> Grid::global_to_local(shared_ptr<Block> in) const
{
  Block r;

  // Size of the local grid, including ghost cells
  unsigned int dxg = info_.dimx_, dyg = info_.dimy_, dzg = info_.dimz_;

  // Size of the local grid, NOT including ghost cells
  unsigned int dx = info_.dimx_no_sd_, dy = info_.dimy_no_sd_, 
    dz = info_.dimz_no_sd_;

  // Starting points of the local grid within the global grid,
  // including ghost cells.
  unsigned int sxg = info_.start_x_, syg = info_.start_y_, 
    szg = info_.start_z_;

  // Starting points of the local grid within the global grid,
  // NOT including ghost cells.
  unsigned int sx = info_.start_x_no_sd_, sy = info_.start_y_no_sd_, 
    sz = info_.start_z_no_sd_;

  r.is_global_ = false;


  // Calculate Block parameters for a new Block in the local grid that
  // DOES NOT include ghost cells.
  if (in->xmin_ >= sx) 
  {
    if (in->xmin_ < sx + dx) {
      r.xmin_ = in->xmin_ - sxg;
      r.faces_[BACK] = true;
    } else {
      r.xmin_ = 0;
      r.faces_[BACK] = false;
      r.has_data_ = false;
    }
  } else {
    r.xmin_ = sx - sxg;
    r.xoffset_ = sx - in->xmin_;
    r.faces_[BACK] = false;
  }

  if (in->ymin_ >= sy) 
  {
    if (in->ymin_ < sy + dy) {
      r.ymin_ = in->ymin_ - syg;
      r.faces_[LEFT] = true;
    } else {
      r.ymin_ = 0;
      r.faces_[LEFT] = false;
      r.has_data_ = false;
    }
  } else {
    r.ymin_ = sy - syg;
    r.yoffset_ = sy - in->ymin_;
    r.faces_[LEFT] = false;
  }

  if (in->zmin_ >= sz)
  {
    if (in->zmin_ < sz + dz) {
      r.zmin_ = in->zmin_ - szg;
      r.faces_[BOTTOM] = true;
    } else {
      r.zmin_ = 0;
      r.faces_[BOTTOM] = false;
      r.has_data_ = false;
    }
  } else {
    r.zmin_ = sz - szg; 
    r.zoffset_ = sz - in->zmin_;
    r.faces_[BOTTOM] = false;
  }

  // Maximums... 
  if (in->xmax_ >= sx) {
    if (in->xmax_ < sx + dx) {
      r.xmax_ = in->xmax_ - sxg;
      r.faces_[FRONT] = true;
    } else {
      r.xmax_ = sx - sxg + dx - 1;
      r.faces_[FRONT] = false;
    }
  } else {
    r.xmax_ = 0;
    r.faces_[FRONT] = false;
  }

  if (in->ymax_ >= sy) {
    if (in->ymax_ < sy + dy) {
      r.ymax_ = in->ymax_ - syg;
      r.faces_[RIGHT] = true;
    } else {
      r.ymax_ = sy - syg + dy - 1;
      r.faces_[RIGHT] = false;
    }
  } else {
    r.ymax_ = 0;
    r.faces_[RIGHT] = false;
  }

  if (in->zmax_ >= sz) {
    if (in->zmax_ < sz + dz) {
      r.zmax_ = in->zmax_ - szg;
      r.faces_[TOP] = true;
    } else {
      r.zmax_ = sz - szg + dz - 1;
      r.faces_[TOP] = false;
    }
  } else {
    r.zmax_ = 0;
    r.faces_[TOP] = false;
  }

  r.xlen_ = r.xmax_ - r.xmin_ + 1;
  r.ylen_ = r.ymax_ - r.ymin_ + 1;
  r.zlen_ = r.zmax_ - r.zmin_ + 1;

  return shared_ptr<Block> (new Block(r));
}

shared_ptr<Block> Grid::global_to_local_ghost(shared_ptr<Block> in) const
{
  Block r;

  // Size of the local grid, including ghost cells
  unsigned int dxg = info_.dimx_, dyg = info_.dimy_, dzg = info_.dimz_;

  // Size of the local grid, NOT including ghost cells
  unsigned int dx = info_.dimx_no_sd_, dy = info_.dimy_no_sd_, 
    dz = info_.dimz_no_sd_;

  // Starting points of the local grid within the global grid,
  // including ghost cells.
  unsigned int sxg = info_.start_x_, syg = info_.start_y_, 
    szg = info_.start_z_;

  // Starting points of the local grid within the global grid,
  // NOT including ghost cells.
  unsigned int sx = info_.start_x_no_sd_, sy = info_.start_y_no_sd_, 
    sz = info_.start_z_no_sd_;

  r.is_global_ = false;


  // Calculate Block parameters for a new Block in the local grid that
  // DOES NOT include ghost cells.
  if (in->xmin_ >= sxg) 
  {
    if (in->xmin_ < sxg + dxg) {
      r.xmin_ = in->xmin_ - sxg;
      r.faces_[BACK] = true;
    } else {
      r.xmin_ = 0;
      r.faces_[BACK] = false;
      r.has_data_ = false;
    }
  } else {
    r.xmin_ = 0;
    r.xoffset_ = sxg - in->xmin_;
    r.faces_[BACK] = false;
  }

  if (in->ymin_ >= syg) 
  {
    if (in->ymin_ < syg + dyg) {
      r.ymin_ = in->ymin_ - syg;
      r.faces_[LEFT] = true;
    } else {
      r.ymin_ = 0;
      r.faces_[LEFT] = false;
      r.has_data_ = false;
    }
  } else {
    r.ymin_ = 0;
    r.yoffset_ = syg - in->ymin_;
    r.faces_[LEFT] = false;
  }

  if (in->zmin_ >= szg)
  {
    if (in->zmin_ < szg + dzg) {
      r.zmin_ = in->zmin_ - szg;
      r.faces_[BOTTOM] = true;
    } else {
      r.zmin_ = 0;
      r.faces_[BOTTOM] = false;
      r.has_data_ = false;
    }
  } else {
    r.zmin_ = 0;
    r.zoffset_ = szg - in->zmin_;
    r.faces_[BOTTOM] = false;
  }

  // Maximums... 
  if (in->xmax_ >= sxg) {
    if (in->xmax_ < sxg + dxg) {
      r.xmax_ = in->xmax_ - sxg;
      r.faces_[FRONT] = true;
    } else {
      r.xmax_ = dxg - 1;
      r.faces_[FRONT] = false;
    }
  } else {
    r.xmax_ = 0;
    r.faces_[FRONT] = false;
  }

  if (in->ymax_ >= syg) {
    if (in->ymax_ < syg + dyg) {
      r.ymax_ = in->ymax_ - syg;
      r.faces_[RIGHT] = true;
    } else {
      r.ymax_ = dyg - 1;
      r.faces_[RIGHT] = false;
    }
  } else {
    r.ymax_ = 0;
    r.faces_[RIGHT] = false;
  }

  if (in->zmax_ >= szg) {
    if (in->zmax_ < szg + dzg) {
      r.zmax_ = in->zmax_ - szg;
      r.faces_[TOP] = true;
    } else {
      r.zmax_ = dzg - 1;
      r.faces_[TOP] = false;
    }
  } else {
    r.zmax_ = 0;
    r.faces_[TOP] = false;
  }

  r.xlen_ = r.xmax_ - r.xmin_ + 1;
  r.ylen_ = r.ymax_ - r.ymin_ + 1;
  r.zlen_ = r.zmax_ - r.zmin_ + 1;

  return shared_ptr<Block> (new Block(r));
}

point Grid::get_size() const
{
  return (*pg_).get_grid_size();
}

point Grid::get_centre() const
{
  return (*pg_).get_grid_centre();
}

/**
 * \bug This doesn't work right if the point is larger than the bounds!!!
 * \bug No way to indicate that the point is actually outside the local grid
 */
grid_point Grid::global_to_local(grid_point p) const
{
  grid_point r;
  
  r.x = (info_.start_x_ > p.x) ? 0 : p.x - info_.start_x_;
  r.y = (info_.start_y_ > p.y) ? 0 : p.y - info_.start_y_;
  r.z = (info_.start_z_ > p.z) ? 0 : p.z - info_.start_z_;

  return r;
}
