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
 */

#ifndef GRID_H
#define GRID_H

#include <mpi.h>

#include <time.h>

#include "Types.hh"
#include "MaterialLib.hh"

class Grid {
 private:
 protected:

  /** 
   * Allocate memory for the grid. Called by setup_grid(). 
   */ 
  virtual void alloc_grid();

  /**
   * Free the memory allocated for the grid. 
   */ 
  virtual void free_grid();

  /**
   * Set up the MPI derived datatypes used to send face information
   * between nodes. 
   */
  void init_datatypes();

  // Global grid size (i.e. all domains)
  unsigned int global_dimx_;
  unsigned int global_dimy_;
  unsigned int global_dimz_;

  // Local grid starting point
  unsigned int start_x_;
  unsigned int start_y_;
  unsigned int start_z_;

  // Local grid size (this sub domain only)
  // If a face must be communicated to a processor, then the
  // dimention will increase by one. If both the top and bottom must
  // be sent to other processors, dimz_ = actual subdomain size +
  // 2. This is overlap which allows us to calculate the interior of
  // the subdomain without having to stop and ask the other
  // processors for data. 
  unsigned int dimx_;
  unsigned int dimy_;
  unsigned int dimz_;

  // Time and space steppings; the distance between each point in the
  // grid.
  delta_t deltax_;
  delta_t deltay_;
  delta_t deltaz_;
  delta_t deltat_;

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

  // A grid is a cube with six faces. Those faces either need to have
  // boundary conditions, or they are subdomain boundaries and they
  // need to be shared with other processors. These arrays tell what
  // to do with each face. 
  //
  // 0 - Front (x = 0, YZ plane)
  // 1 - Back (x = dimx, YZ plane)
  // 2 - Left (y = 0, XZ plane)
  // 3 - Right (y = dimy, XZ plane)
  // 4 - Bottom (z = 0, XY plane)
  // 5 - Top (z = dimz, XY plane)
  //
  BoundaryCondition face_bc_[6]; // Boundary condition to apply
  int face_rank_[6]; // Rank of processor to talk to about this
		     // interface. 

  // Derived MPI data types for sending data around. 
  MPI_Datatype xy_plane_;
  MPI_Datatype yz_plane_;
  MPI_Datatype xz_plane_;

  // Accessing Z is contiguous, but accessing X and Y coordinates
  // requires the use of strided vectors. 
  MPI_Datatype x_vector_;
  MPI_Datatype y_vector_;
  MPI_Datatype z_vector_;

 public:
  Grid();
  virtual ~Grid();

  /**
   * Compute the next time step of the fields. The update equations
   * used here are straight out of Taflove. 
   */
  virtual void update_fields();

  /**
   * Apply the boundary conditions to the faces
   */
  virtual void apply_boundaries();

  /**
   * Calculate the material constants from the given material library
   *
   * @param matlib the material library to load the materials from 
   */
  virtual void load_materials(MaterialLib &matlib);

  /**
   * Deallocate the memory used to store material coeffcients and so on. 
   */
  virtual void free_material();

  /**
   * Set up the, define the size of the global grid, and the size of
   * the subdomain this grid actually represents. 
   *
   * @param global_x number of cells in the global grid along the x dimension
   * @param global_y number of cells in the global grid along the y dimension
   * @param global_z number of cells in the global grid along the z dimension
   *
   * @param start_x the starting cell of the subdomain along x
   * @param start_y the starting cell of the subdomain along y
   * @param start_z the starting cell of the subdomain along z
   *
   * @param x the size of the subdomain in the x dimension
   * @param y the size of the subdomain in the y dimension
   * @param z the size of the subdomain in the z dimension
   *
   * @param deltax the spacing between nodes along x
   * @param deltay the spacing between nodes along y
   * @param deltaz the spacing between nodes along z
   * @param deltat the time between time steps
   */ 
  virtual void setup_grid(unsigned int global_x, unsigned int global_y, 
                          unsigned int global_z, 
                          unsigned int start_x, unsigned int start_y, 
                          unsigned int start_z, 
                          unsigned int x, unsigned int y, unsigned int z, 
                          delta_t deltax, delta_t deltay, delta_t deltaz,
                          delta_t deltat);

  /**
   * Define geometry in the grid (i.e. assign material indicies to
   * grid points). All definitions are done in global coordinates. It
   * is the grid's job to calculate the region within the subdomain
   * and assign the material indicies appropriatly. 
   *
   * @param x_start The starting cell of the box in x
   * @param x_end The ending cell of the box in x
   * @param y_start The starting cell of the box in y
   * @param y_end The ending cell of the box in y
   * @param z_start The starting cell of the box in z
   * @param z_end The ending cell of the box in z
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
   * Set the boundary condition on one of the faces of this grid. 
   *
   * @param face the face to assign the boundary to. One of FRONT, BACK,
   * LEFT, RIGHT, BOTTOM, TOP as defined in Types.hh
   *
   * @param bc the boundary condition to apply, one of SUBDOMAIN
   * (meaning that information is exchanged with another node at this
   * face), EWALL, HWALL, ESYM, HSYM, IMPEDANCE, or PML, as defined in
   * Types.hh
   */ 
  void set_boundary(unsigned int face, BoundaryCondition bc);

  /**
   * Set the rank of the processor this face needs to be shared with. 
   *
   * @param face the face to assign the node rank to. One of FRONT, BACK,
   * LEFT, RIGHT, BOTTOM, TOP as defined in Types.hh
   *
   * @param rank the rank of the node to share this face with. 
   */ 
  void set_face_rank(unsigned int face, int rank);

  // Accessors (these should be inline, but then I would have to
  // rearrange definitions, and I'm too lazy right now. 
  
  /**
   * Returns the global size of the x dimension.
   *
   * @return global grid x size
   */
  unsigned int get_gdx();
  /**
   * Returns the global size of the y dimension.
   *
   * @return global grid y size
   */
  unsigned int get_gdy();
  /**
   * Returns the global size of the z dimension.
   *
   * @return global grid z size
   */
  unsigned int get_gdz();

  /**
   * Returns the local start of the x dimension.
   *
   * @return local grid x start
   */
  unsigned int get_lsx();
  /**
   * Returns the local start of the y dimension.
   *
   * @return local grid y start
   */
  unsigned int get_lsy();
  /**
   * Returns the local start of the z dimension.
   *
   * @return local grid z start
   */
  unsigned int get_lsz();


  /**
   * Returns the local size of the x dimension.
   *
   * @return local grid x size
   */
  unsigned int get_ldx();
  /**
   * Returns the local size of the y dimension.
   *
   * @return local grid y size
   */
  unsigned int get_ldy();
  /**
   * Returns the local size of the z dimension.
   *
   * @return local grid z size
   */
  unsigned int get_ldz();


  /**
   * Returns the spacing in the x direction (mm)
   *
   * @return the x spacing in mm
   */
  delta_t get_deltax();
  /**
   * Returns the spacing in the y direction (mm)
   *
   * @return the y spacing in mm
   */
  delta_t get_deltay();
  /**
   * Returns the spacing in the z direction (mm)
   *
   * @return the z spacing in mm
   */
  delta_t get_deltaz();
  /**
   * Returns the spacing in time (s)
   *
   * @return the time spacing in seconds
   */
  delta_t get_deltat();

  // Assign and get values of field components

  /**
   * Assign a value to Ex at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @param val the value to set ex to
   */
  void set_ex(unsigned int x, unsigned int y, 
                     unsigned int z, field_t val);

  /**
   * Assign a value to Ey at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @param val the value to set ey to
   */
  void set_ey(unsigned int x, unsigned int y, 
                     unsigned int z, field_t val);

  /**
   * Assign a value to Ez at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @param val the value to set ez to
   */
  void set_ez(unsigned int x, unsigned int y, 
                     unsigned int z, field_t val);

  /**
   * Assign a value to Hx at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @param val the value to set hx to
   */
  void set_hx(unsigned int x, unsigned int y, 
                     unsigned int z, field_t val);

  /**
   * Assign a value to Hy at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @param val the value to set hy to
   */
  void set_hy(unsigned int x, unsigned int y, 
                     unsigned int z, field_t val);

  /**
   * Assign a value to Hz at some point in space.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param z the z coordinate
   * @param val the value to set hz to
   */
  void set_hz(unsigned int x, unsigned int y, 
                     unsigned int z, field_t val);
};

#endif // GRID_H

