/* 
   phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004 Matt Hughes <mhughe@uvic.ca>

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

/** \class Grid
 * \brief Base class for the grid. 
 * 
 * Data for the grid is defined here.  The basic update equations are
 * defined by this class, but they can be subclassed to take advantage
 * of certian things only available on certian arches (SSE2, AltiVec),
 * or to provide different types of materials, such as Debye, Lorentz,
 * etc.
 *
 * Since this is an MPI program, and each processor only deals with
 * part of domain. This Grid class is owned by the processor, and we
 * assume that there is only one grid per processor. If that's not the
 * case, well, something is probably broken.
 *
 * \bug Since we have this define mode now, and a lot of the functions
 * have to check that flag, there should be protected virtual helper
 * functions that can be overridden to do the actual work, while the
 * public functions are nonvirtual and just perform the define mode
 * check and call the protected helper functions. Saves derived
 * classes from having to duplicate code or messing up define mode. 
 */

#ifndef GRID_H
#define GRID_H

#include "config.h" // must come before assert.h

#include <mpi.h>

#include <time.h>
#include <assert.h>

#include "Types.hh"
#include "MaterialLib.hh"
#include "GridInfo.hh"
#include "Geometry.hh"

#include <map>

using namespace std;

class Grid {
  friend class UPml;
  friend class UPmlCommon;
  friend class Pml; // So that PML update equations can access the
                    // field pointers quickly; and also so that they
                    // can make a giant mess of them. 
  friend class PmlCommon;

 private:
 protected:

  /**
   * Free the memory allocated for the grid. 
   */ 
  virtual void free_grid();

  /**
   * Set up the MPI derived datatypes used to send face information
   * between nodes. 
   */
  void init_datatypes();

  /**
   * Free the MPI derived datatypes
   */
  void free_datatypes();

  // Grid size information
  GridInfo info_;

  // These help simply things when dealing with non-zero thickness
  // boundary conditions (PML), and help simplify update equations by
  // removing the need to do math with regard to what can be
  // updated. Set up in set_define_mode(false).

  region_t update_ex_r_; /**< Region over which ex update is applied. */
  region_t update_ey_r_; /**< Region over which ey update is applied. */
  region_t update_ez_r_; /**< Region over which ez update is applied. */

  region_t update_hx_r_; /**< Region over which hx update is applied. */
  region_t update_hy_r_; /**< Region over which hy update is applied. */
  region_t update_hz_r_; /**< Region over which hz update is applied. */

  // Number of materials we know about (0 is PEC)
  unsigned int num_materials_;

  // E Field Material Coefficients
  mat_coef_t *Ca_;
  mat_coef_t *Cbx_;
  mat_coef_t *Cby_;
  mat_coef_t *Cbz_;
  
  // H Field Coefficients
  mat_coef_t *Da_;
  mat_coef_t *Dbx_;
  mat_coef_t *Dby_;
  mat_coef_t *Dbz_;

  // Material library
  MaterialLib material_lib_;

  /**
   * Field data: one big vector, where a 3d coordinate is indexed
   * using idx = z + (y + x*dimy) * dimz. This is done this way so
   * that is is easy to use MPI derived data types to move data
   * around. It should also allow for a increase in speed, since data
   * access requires only one pointer dereference. 
   */
  field_t *ex_;
  field_t *ey_;
  field_t *ez_;
  field_t *hx_;
  field_t *hy_;
  field_t *hz_;

  // The material for each point in the grid. This is an index into
  // the material arrays, Ca, Cbx, etc. 
  unsigned int *material_;

  // Derived MPI data types for sending data around. 
  MPI_Datatype xy_plane_;
  MPI_Datatype yz_plane_;
  MPI_Datatype xz_plane_;

  // Accessing Z is contiguous, but accessing X and Y coordinates
  // requires the use of strided vectors. 
  MPI_Datatype x_vector_;
  MPI_Datatype y_vector_;
  MPI_Datatype z_vector_;

  bool types_alloced_; 

  /**
   * This is true when the grid is in define mode. In define mode,
   * the grid size can be changed, material definitions can be
   * changed, boundary conditions can be set, etc. The update
   * equations cannot be run in define mode. 
   */
  bool define_;

  /**
   * An array of pointers to geometry objects. Since this will be used
   * inside the update equation loops, it is an array rather than a
   * vector to avoid iterator construction overhead and so
   * on. Hopefully the virtual function calls won't be too much of a
   * problem. 
   */ 
  Geometry **geometries_;
  unsigned int num_geoms_;

  /** 
   * A sancturary for boundary condition common data. The valid int's
   * that can be used as keys are enumerated in Types.hh.in.
   */
  map<GridAuxData, void *> auxdata_; 

  /**
   * Store some auxiliary data. 
   */
  inline void add_auxdata(GridAuxData n, void *ptr)
  {
    auxdata_[n] = ptr;
  }

  /**
   * Retrieves some auxiliry data. Returns null if the requested data
   * is not available.
   */
  inline void *get_auxdata(GridAuxData n)
  {
    map<GridAuxData, void *>::iterator iter = auxdata_.find(n);
    if (iter == auxdata_.end())
      return 0;
    else
      return iter->second;
  }

  /**
   * Compute the update equatations for the Ex field component. 
   *
   * @param update_r the region to update. A param in case a boundary
   * condition wants to call it. 
   */
  virtual void update_ex(region_t update_r);

  /**
   * Compute the update equatations for the Ey field component. 
   *
   * @param update_r the region to update. A param in case a boundary
   * condition wants to call it. 
   */
  virtual void update_ey(region_t update_r);

  /**
   * Compute the update equatations for the Ez field component. 
   *
   * @param update_r the region to update. A param in case a boundary
   * condition wants to call it. 
   */
  virtual void update_ez(region_t update_r);

  /**
   * Compute the update equatations for the Hx field component. 
   *
   * @param update_r the region to update. A param in case a boundary
   * condition wants to call it. 
   */
  virtual void update_hx(region_t update_r);

  /**
   * Compute the update equatations for the Hy field component. 
   *
   * @param update_r the region to update. A param in case a boundary
   * condition wants to call it. 
   */
  virtual void update_hy(region_t update_r);

  /**
   * Compute the update equatations for the Hz field component. 
   *
   * @param update_r the region to update. A param in case a boundary
   * condition wants to call it. 
   */
  virtual void update_hz(region_t update_r);

  /** 
   * Allocate memory for the grid. This function is called when
   * set_define_mode() turns define mode off.  
   */ 
  virtual void alloc_grid();

  /**
   * This is called by set_define_mode(true), and it should add any
   * data that needs to be exchanged across subdomain boundaries to 
   * the SubdomainBc object. 
   *
   * @param sd the subdomain boundary to add data to be exchanged to.
   * @param face the face this subdomain boundary is on. 
   */
  virtual void setup_subdomain_data(SubdomainBc *sd, Face face);

 public:
  Grid();
  virtual ~Grid();

  /**
   * Point Index: Calculate the index in the arrays of a 3d
   * coordinate. ALWAYS USE THIS FUNCTION, in case I change the way
   * things are organized for some reason. It's inline, so it should
   * compile out.
   *
   * @param x
   * @param y
   * @param z
   * @param an index into the field component and material arrays. 
   */
  inline unsigned int pi(unsigned int x, unsigned int y, 
                         unsigned int z) const
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    return z + (y + x*info_.dimy_) * info_.dimz_;
  }

  /**
   * Copy constructor. Pretty much everything is copied, but the grid
   * is condidered uninitialized. No memory is allocated, grid data
   * pointers are set to zero. 
   */
  Grid(const Grid &rhs);

  /**
   * Assignment operator. Pretty much everything is copied, but the
   * grid is condidered uninitialized. No memory is allocated, grid
   * data pointers are set to zero.
   */
  const Grid &operator=(const Grid &rhs);

  /**
   * Return whether the grid is in define mode or not. 
   */
  inline bool get_define_mode() const
  {
    return define_;
  }

  /**
   * Turn define mode on or off. If you are turning off define mode,
   * a number of sanity checks and stability checks are made to
   * ensure that the settings make sense and can be solved. 
   *
   * @param d a boolean, false to turn off define mode, true to turn
   * it on. 
   */
  virtual void set_define_mode(bool d);

  /**
   * Compute the next time step of electric field. This is a convenience
   * function which calls the individual update functions, and can
   * only be called when the grid is not in define mode. 
   */
  void update_e_field();

  /**
   * Compute the next time step of magnetic field. This is a convenience
   * function which calls the individual update functions, and can
   * only be called when the grid is not in define mode. 
   */
  void update_h_field();

  /**
   * Apply the boundary conditions to the faces. This function only
   * has an effect when the grid is not in define mode. 
   */
  virtual void apply_boundaries(FieldType type);

  /**
   * Calculate the material constants from the given material
   * library. This function can only be used in define mode. 
   *
   * @param matlib the material library to load the materials from 
   */
  virtual void load_materials(MaterialLib &matlib);

  /**
   * Store references to the geometry objects so that more specialized
   * dispersions can store data in them.
   */ 
  virtual void load_geometries(vector<Geometry *> &geoms);

  /**
   * Deallocate the memory used to store material coeffcients and so
   * on. This function can only be used in define mode. 
   */
  virtual void free_material();

  /**
   * Set up the, define the size of the global grid, and the size of
   * the subdomain this grid actually represents. This function can
   * only be used in define mode. 
   *
   * @param info the GridInfo object containing the global and local
   * grid sizes, cell sizes, time step size, etc. 
   */ 
  virtual void setup_grid(const GridInfo &info);

  /**
   * Retrieve information about the grid's size, space and time step
   * size, and boundary conditions. 
   *
   * @return a reference to the GridInfo object 
   */
  inline GridInfo& get_grid_info() 
  {
    return info_;
  }

  /** 
   * Returns a reference to the boundary condition object for a face. 
   *
   * @return ref to a BoundaryCond
   */
  inline const BoundaryCond& get_boundary(Face face) const
  {
    return const_cast<const BoundaryCond &>(info_.get_boundary(face));
  }

  /**
   * Convert global coordinate to local grid coordinates. Generally
   * excitations and geometry specifications are given in terms of
   * the global FDTD grid. Since each processor only knows about part
   * of that grid, it is necessary to convert global coordinates to
   * local ones that can be used. 
   *
   * @param r the global region to convert.
   * @param no_ol true if the subdomain overlap should NOT be included (for Results)
   *
   * @return region in local coordinate. 
   */
  region_t global_to_local(region_t r, bool no_ol = false) const;

  /**
   * Convert a global point to local coordinates.
   *
   * @param p the point to convert
   * @return the point in global coordinates
   */
  point_t global_to_local(point_t p) const;

  /**
   * Convert global coordinate to local grid coordinates. Generally
   * excitations and geometry specifications are given in terms of
   * the global FDTD grid. Since each processor only knows about part
   * of that grid, it is necessary to convert global coordinates to
   * local ones that can be used. 
   *
   * @param x_start The starting cell of the box in x
   * @param x_stop The ending cell of the box in x
   * @param y_start The starting cell of the box in y
   * @param y_stop The ending cell of the box in y
   * @param z_start The starting cell of the box in z
   * @param z_stop The ending cell of the box in z
   *
   * @return region_t in local coordinate. 
   */
  region_t global_to_local(unsigned int x_start, unsigned int x_stop, 
                           unsigned int y_start, unsigned int y_stop, 
                           unsigned int z_start, unsigned int z_stop,
                           bool no_ol = false) const;

  /**
   * Define geometry in the grid (i.e. assign material indicies to
   * grid points). All definitions are done in global coordinates. It
   * is the grid's job to calculate the region within the subdomain
   * and assign the material indicies appropriatly. This function can
   * only be used in define mode. 
   *
   * @param x_start The starting cell of the box in x
   * @param x_stop The ending cell of the box in x
   * @param y_start The starting cell of the box in y
   * @param y_stop The ending cell of the box in y
   * @param z_start The starting cell of the box in z
   * @param z_stop The ending cell of the box in z
   *
   * @param mat_index The material index to use in this region. 0 is
   * perfect electric conductor, 1 and up are ordered as in the
   * material library.
   */
  //virtual void define_box(unsigned int x_start, unsigned int x_stop, 
  //                        unsigned int y_start, unsigned int y_stop, 
  //                        unsigned int z_start, unsigned int z_stop, 
  //                        unsigned int mat_index);

  /**
   * Returns the global size of the x dimension.
   *
   * @return global grid x size
   */
  inline unsigned int get_gdx() const
  {
    return info_.global_dimx_;
  }

  /**
   * Returns the global size of the y dimension.
   *
   * @return global grid y size
   */
  inline unsigned int get_gdy() const
  {
    return info_.global_dimy_;
  }

  /**
   * Returns the global size of the z dimension.
   *
   * @return global grid z size
   */
  inline unsigned int get_gdz() const
  {
    return info_.global_dimz_;
  }

  /**
   * Returns the local start of the x dimension.
   *
   * @return local grid x start
   */
  inline unsigned int get_lsx() const
  {
    return info_.start_x_; // no_sd_;
  }

  /**
   * Returns the local start of the y dimension.
   *
   * @return local grid y start
   */
  inline unsigned int get_lsy() const
  {
    return info_.start_y_; // no_sd_;
  }

  /**
   * Returns the local start of the z dimension.
   *
   * @return local grid z start
   */
  inline unsigned int get_lsz() const
  {
    return info_.start_z_; // no_sd_;
  }



  /**
   * Returns the local size of the x dimension.
   *
   * @return local grid x size
   */
  inline unsigned int get_ldx() const 
  {
    return info_.dimx_no_sd_;
  }

  /**
   * Returns the local size of the y dimension.
   *
   * @return local grid y size
   */
  inline unsigned int get_ldy() const 
  {
    return info_.dimy_no_sd_;
  }

  /**
   * Returns the local size of the z dimension.
   *
   * @return local grid z size
   */
  inline unsigned int get_ldz() const
  {
    return info_.dimz_no_sd_;
  }

  /**
   * Returns the local size of the x dimension, including any
   * subdomain overlap, if it exists.
   *
   * @return local grid x size
   */
  inline unsigned int get_ldx_sd() const 
  {
    return info_.dimx_;
  }

  /**
   * Returns the local size of the y dimension, including any
   * subdomain overlap, if it exists.
   *
   * @return local grid y size
   */
  inline unsigned int get_ldy_sd() const 
  {
    return info_.dimy_;
  }

  /**
   * Returns the local size of the z dimension, including any
   * subdomain overlap, if it exists.
   *
   * @return local grid z size
   */
  inline unsigned int get_ldz_sd() const
  {
    return info_.dimz_;
  }

  /**
   * Returns the spacing in the x direction (mm)
   *
   * @return the x spacing in mm
   */
  inline delta_t get_deltax() const 
  {
    return info_.deltax_;
  }

  /**
   * Returns the spacing in the y direction (mm)
   *
   * @return the y spacing in mm
   */
  inline delta_t get_deltay() const 
  {
    return info_.deltay_;
  }

  /**
   * Returns the spacing in the z direction (mm)
   *
   * @return the z spacing in mm
   */
  inline delta_t get_deltaz() const 
  {
    return info_.deltaz_;
  }

  /**
   * Returns the spacing in time (s)
   *
   * @return the time spacing in seconds
   */
  inline delta_t get_deltat() const 
  {
    return info_.deltat_;
  }

  // Assign and get values of field components

  /**
   * Assign a value to Ex at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @param val the value to set ex to
   */
  inline void set_ex(unsigned int x, unsigned int y, 
                     unsigned int z, field_t val)
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    ex_[pi(x, y, z)] = val;
  }

  /**
   * Assign a value to Ey at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @param val the value to set ey to
   */
  inline void set_ey(unsigned int x, unsigned int y, 
                     unsigned int z, field_t val)
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    ey_[pi(x, y, z)] = val;
  }


  /**
   * Assign a value to Ez at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @param val the value to set ez to
   */
  inline void set_ez(unsigned int x, unsigned int y, 
                     unsigned int z, field_t val)
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    ez_[pi(x, y, z)] = val;
  }


  /**
   * Assign a value to Hx at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @param val the value to set hx to
   */
  inline void set_hx(unsigned int x, unsigned int y, 
                     unsigned int z, field_t val)
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    hx_[pi(x, y, z)] = val;
  }

  /**
   * Assign a value to Hy at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @param val the value to set hy to
   */
  inline void set_hy(unsigned int x, unsigned int y, 
                     unsigned int z, field_t val)
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    hy_[pi(x, y, z)] = val;
  }

  /**
   * Assign a value to Hz at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @param val the value to set hz to
   */
  inline void set_hz(unsigned int x, unsigned int y, 
                     unsigned int z, field_t val)
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    hz_[pi(x, y, z)] = val;
  }


  /**
   * Return Ex at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @return the value of ex
   */
  inline field_t get_ex(unsigned int x, unsigned int y, 
                        unsigned int z) const
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    return ex_[pi(x, y, z)];
  }

  /**
   * Return Ey at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @return the value of ey
   */
  inline field_t get_ey(unsigned int x, unsigned int y, 
                        unsigned int z) const
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    return ey_[pi(x, y, z)];
  }


  /**
   * Return Ez at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @return the value of ez
   */
  inline field_t get_ez(unsigned int x, unsigned int y, 
                        unsigned int z) const
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    return ez_[pi(x, y, z)];
  }


  /**
   * Return Hx at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @return the value of hx
   */
  inline field_t get_hx(unsigned int x, unsigned int y, 
                        unsigned int z) const
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    return hx_[pi(x, y, z)];
  }

  /**
   * Return Hy at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @return the value of hy
   */
  inline field_t get_hy(unsigned int x, unsigned int y, 
                        unsigned int z) const
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    return hy_[pi(x, y, z)];
  }

  /**
   * Return Hz at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @return the value of hz
   */
  inline field_t get_hz(unsigned int x, unsigned int y, 
                        unsigned int z) const
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    return hz_[pi(x, y, z)];
  }

  /**
   * Return the MPI derived data type for the xy plane. 
   *
   * @param face the face to return the data type for
   * @return MPI_Datatype
   */
  MPI_Datatype get_plane_dt(Face face) const;

  /**
   * Return the MPI Derived data type for an X vector
   * @return MPI_Datatype
   */
  inline MPI_Datatype get_x_vector_dt() const
  {
    return x_vector_;
  }

  /**
   * Return the MPI Derived data type for an Y vector
   * @return MPI_Datatype
   */
  inline MPI_Datatype get_y_vector_dt() const
  {
    return y_vector_;
  }

  /**
   * Return the MPI Derived data type for an Z vector
   * @return MPI_Datatype
   */
  inline MPI_Datatype get_z_vector_dt() const
  {
    return z_vector_;
  }

  /**
   * Return a pointer to the start of a face. DANGER!! Clients must
   * take care when using this pointer! 
   *
   * This is intended to be used by SubdomainBc for sending and
   * recieving data between ranks where the grids overlap. 
   *
   * @param face The face of interest
   * @param comp The component of interest
   * @param offset An offset from the face, defaults to zero. Must be
   * greater than zero; results in a pointer to a face inside the
   * grid. 
   * @return a pointer to the field component at the specified face
   */
  field_t *get_face_start(Face face, FieldComponent comp,
                          unsigned int offset = 0) const;

  /**
   * Return a pointer to the start of a face. DANGER!! Clients must
   * take care when using this pointer! 
   *
   * This is intended for use by PlaneResult for getting results out. 
   *
   * @param face The face of interest
   * @param comp The component of interest
   * @param p A point indicating where the plane intersects the grid
   * @return a pointer to the field component at the specified face
   */
  field_t *get_face_start(Face face, FieldComponent comp,
                          point_t p) const;

  /**
   * Return a pointer to a point in the grid for a particular field
   * component. DANGER!! Clients must take care when using this
   * pointer!
   *
   * @param point the point in the grid to get the pointer at
   * @param field_comp the field component to return the pointer for
   * @return a pointer, which may be null if something went wrong. 
   */
  const field_t *get_pointer(point_t point, 
                             FieldComponent field_comp) const;

  /**
   * Assign a material index a some point in space
   *
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @param mid the material id to assign
   */
  inline void set_material(unsigned int x, unsigned int y, 
                           unsigned int z, unsigned int mid)
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    material_[pi(x, y, z)] = mid;
  }

  /**
   * Return a material index a some point in space
   *
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @return the material id
   */
  inline unsigned int get_material(unsigned int x, unsigned int y, 
                                   unsigned int z)
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    return material_[pi(x, y, z)];
  }

  /**
   * Return the material library
   */ 
  inline const MaterialLib &get_material_lib() const
  {
    return material_lib_;
  }

  /**
   * Returns the number of nodes processed by the update equations. 
   */ 
  unsigned int get_num_updated_nodes()
  {
    return info_.dimx_ * info_.dimy_ * info_.dimz_; 
      //  update_ex_r_.xmax - update_ex_r_.xmin + 
//        update_ey_r_.xmax - update_ey_r_.xmin + 
//        update_ez_r_.xmax - update_ez_r_.xmin + 
//        update_hx_r_.xmax - update_hx_r_.xmin + 
//        update_hy_r_.xmax - update_hy_r_.xmin + 
//        update_hz_r_.xmax - update_hz_r_.xmin;
  }

};

#endif // GRID_H

