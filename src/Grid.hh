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
 * \bug A copy constructor, assignment operator, and the other thing
 * are needed so that pointers aren't inadvertantly copies around by
 * accident. 
 *
 * \bug Sanity checks on stability are required.
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

class Grid {
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

  // Grid size information
  GridInfo info_;

  // The region over which the Grid update equations are
  // applied. Just the size given in info_ minus the offsets. This is
  // calculated when we leave define mode. 
  region_t update_r_;

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

  // Field data
  field_t ***ex_;
  field_t ***ey_;
  field_t ***ez_;
  field_t ***hx_;
  field_t ***hy_;
  field_t ***hz_;

  // The material for each point in the grid. This is an index into
  // the material arrays, Ca, Cbx, etc. 
  unsigned int ***material_;

  // Derived MPI data types for sending data around. 
  MPI_Datatype xy_plane_;
  MPI_Datatype yz_plane_;
  MPI_Datatype xz_plane_;

  // Accessing Z is contiguous, but accessing X and Y coordinates
  // requires the use of strided vectors. 
  MPI_Datatype x_vector_;
  MPI_Datatype y_vector_;
  MPI_Datatype z_vector_;

  /**
   * This is true when the grid is in define mode. In define mode,
   * the grid size can be changed, material definitions can be
   * changed, boundary conditions can be set, etc. The update
   * equations cannot be run in define mode. 
   */
  bool define_;


  /**
   * Compute the update equatations for the Ex field component. 
   */
  virtual void update_ex();

  /**
   * Compute the update equatations for the Ey field component. 
   */
  virtual void update_ey();

  /**
   * Compute the update equatations for the Ez field component. 
   */
  virtual void update_ez();

  /**
   * Compute the update equatations for the Hx field component. 
   */
  virtual void update_hx();

  /**
   * Compute the update equatations for the Hy field component. 
   */
  virtual void update_hy();

  /**
   * Compute the update equatations for the Hz field component. 
   */
  virtual void update_hz();

 public:
  Grid();
  virtual ~Grid();

  /**
   * Turn define mode on or off. If you are turning off define mode,
   * a number of sanity checks and stability checks are made to
   * ensure that the settings make sense and can be solved. 
   *
   * @param d a boolean, false to turn off define mode, true to turn
   * it on. 
   */
  void set_define_mode(bool d);

  /**
   * Compute the next time step of the fields. This is a convenience
   * function which calls the individual update functions, and can
   * only be called when the grid is not in define mode. 
   */
  void update_fields();

  /**
   * Apply the boundary conditions to the faces. This function only
   * has an effect when the grid is not in define mode. 
   */
  virtual void apply_boundaries();

  /**
   * Calculate the material constants from the given material
   * library. This function can only be used in define mode. 
   *
   * @param matlib the material library to load the materials from 
   */
  virtual void load_materials(MaterialLib &matlib);

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
   * @return a copy of the GridInfo object 
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
  inline BoundaryCond& get_boundary(Face face)
  {
    return info_.get_boundary(face);
  }

  /** 
   * Allocate memory for the grid. This function can only be used in
   * define mode. 
   */ 
  virtual void alloc_grid();

  /**
   * Convert global coordinate to local grid coordinates. Generally
   * excitations and geometry specifications are given in terms of
   * the global FDTD grid. Since each processor only knows about part
   * of that grid, it is necessary to convert global coordinates to
   * local ones that can be used. 
   *
   * @param r the global region to convert.
   * @return region in local coordinate. 
   */
  region_t global_to_local(region_t r);

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
                         unsigned int z_start, unsigned int z_stop);

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
  void define_box(unsigned int x_start, unsigned int x_stop, 
                  unsigned int y_start, unsigned int y_stop, 
                  unsigned int z_start, unsigned int z_stop, 
                  unsigned int mat_index);

  /**
   * Returns the global size of the x dimension.
   *
   * @return global grid x size
   */
  inline unsigned int get_gdx()
  {
    return info_.global_dimx_;
  }

  /**
   * Returns the global size of the y dimension.
   *
   * @return global grid y size
   */
  inline unsigned int get_gdy()
  {
    return info_.global_dimy_;
  }

  /**
   * Returns the global size of the z dimension.
   *
   * @return global grid z size
   */
  inline unsigned int get_gdz()
  {
    return info_.global_dimz_;
  }

  /**
   * Returns the local start of the x dimension.
   *
   * @return local grid x start
   */
  inline unsigned int get_lsx()
  {
    return info_.start_x_;
  }

  /**
   * Returns the local start of the y dimension.
   *
   * @return local grid y start
   */
  inline unsigned int get_lsy()
  {
    return info_.start_y_;
  }

  /**
   * Returns the local start of the z dimension.
   *
   * @return local grid z start
   */
  inline unsigned int get_lsz()
  {
    return info_.start_z_;
  }



  /**
   * Returns the local size of the x dimension.
   *
   * @return local grid x size
   */
  inline unsigned int get_ldx()
  {
    return info_.dimx_;
  }

  /**
   * Returns the local size of the y dimension.
   *
   * @return local grid y size
   */
  inline unsigned int get_ldy()
  {
    return info_.dimy_;
  }

  /**
   * Returns the local size of the z dimension.
   *
   * @return local grid z size
   */
  inline unsigned int get_ldz()
  {
    return info_.dimz_;
  }


  /**
   * Returns the spacing in the x direction (mm)
   *
   * @return the x spacing in mm
   */
  inline delta_t get_deltax()
  {
    return info_.deltax_;
  }

  /**
   * Returns the spacing in the y direction (mm)
   *
   * @return the y spacing in mm
   */
  inline delta_t get_deltay()
  {
    return info_.deltay_;
  }

  /**
   * Returns the spacing in the z direction (mm)
   *
   * @return the z spacing in mm
   */
  inline delta_t get_deltaz()
  {
    return info_.deltaz_;
  }

  /**
   * Returns the spacing in time (s)
   *
   * @return the time spacing in seconds
   */
  inline delta_t get_deltat()
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
    ex_[x][y][z] = val;
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
    ey_[x][y][z] = val;
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
    ez_[x][y][z] = val;
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
    hx_[x][y][z] = val;
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
    hy_[x][y][z] = val;
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
    hz_[x][y][z] = val;
  }


  /**
   * Return Ex at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @return the value of ex
   */
  inline field_t get_ex(unsigned int x, unsigned int y, 
                        unsigned int z)
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    return ex_[x][y][z];
  }

  /**
   * Return Ey at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @return the value of ey
   */
  inline field_t get_ey(unsigned int x, unsigned int y, 
                        unsigned int z)
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    return ey_[x][y][z];
  }


  /**
   * Return Ez at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @return the value of ez
   */
  inline field_t get_ez(unsigned int x, unsigned int y, 
                        unsigned int z)
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    return ez_[x][y][z];
  }


  /**
   * Return Hx at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @return the value of hx
   */
  inline field_t get_hx(unsigned int x, unsigned int y, 
                        unsigned int z)
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    return hx_[x][y][z];
  }

  /**
   * Return Hy at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @return the value of hy
   */
  inline field_t get_hy(unsigned int x, unsigned int y, 
                        unsigned int z)
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    return hy_[x][y][z];
  }

  /**
   * Return Hz at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @return the value of hz
   */
  inline field_t get_hz(unsigned int x, unsigned int y, 
                        unsigned int z)
  {
    assert(x < info_.dimx_ && y < info_.dimy_ && z < info_.dimz_);
    return hz_[x][y][z];
  }

};

#endif // GRID_H

