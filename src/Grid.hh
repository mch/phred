// Base class for the grid. Data for the grid is defined here.
// The basic update equations are defined by this class, but
// they can be subclassed to take advantage of certian things
// only available on certian arches (SSE2, AltiVec), or to 
// provide different types of materials, such as Debye, Lorentz, 
// etc.

// Since this is an MPI program, and each processor only deals with
// part of domain. This Grid class is owned by the processor, and we
// assume that there is only one grid per processor. If that's not
// the case, well, something is probably broken. 

#ifndef GRID
#define GRID

#include <mpi.h>

#include <time.h>

#include "Types.hh"
#include "MaterialLib.hh"

class Grid {
 private:
 protected:
  void alloc_grid();
  void free_grid();
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

  // Field data, {x,y,z}
  field_t ***ex_;
  field_t ***ey_;
  field_t ***ez_;
  field_t ***hx_;
  field_t ***hy_;
  field_t ***hz_;

  // Running sums of field data, for dispersive material. These aren't
  // allocated if they aren't needed. Most materials of interest are
  // non-magnetic, so we won't bother with the H components for now. 
  field_t ***ex_sum_;
  field_t ***ey_sum_;
  field_t ***ez_sum_;

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

  // These can be overridden; perhaps for AltiVec or SSE2, or for
  // something simple and fast. 
  virtual void update_e();
  virtual void update_h();

  // Calculate the material constants from the given material library
  void load_materials(MaterialLib &matlib);
  void free_material();

  // Grid actions
  void setup_grid(unsigned int global_x, unsigned int global_y, 
                  unsigned int global_z, 
                  unsigned int start_x, unsigned int start_y, 
                  unsigned int start_z, 
                  unsigned int x, unsigned int y, unsigned int z, 
                  delta_t deltax, delta_t deltay, delta_t deltaz,
                  delta_t deltat);

  // Define geometry in the grid (i.e. assign material indicies to
  // grid points). All definitions are done in global coordinates. It
  // is the grid's job to calculate the region within the subdomain
  // and assign the material indicies appropriatly. 
  void define_box(unsigned int x_start, unsigned int x_stop, 
                  unsigned int y_start, unsigned int y_stop, 
                  unsigned int z_start, unsigned int z_stop, 
                  unsigned int mat_index);

  // Boundary condition and subdomain boundary operations. 
  void set_boundary(unsigned int face, BoundaryCondition bc);

  // Set the rank of the processor this face needs to be shared with. 
  void set_face_rank(unsigned int face, int rank);

  // Accessors
  
  /**
   * Returns the global size of the x dimension.
   */
  unsigned int get_gdx();
  /**
   * Returns the global size of the y dimension.
   */
  unsigned int get_gdy();
  /**
   * Returns the global size of the z dimension.
   */
  unsigned int get_gdz();

  /**
   * Returns the local start of the x dimension.
   */
  unsigned int get_lsx();
  /**
   * Returns the local start of the y dimension.
   */
  unsigned int get_lsy();
  /**
   * Returns the local start of the z dimension.
   */
  unsigned int get_lsz();


  /**
   * Returns the local size of the x dimension.
   */
  unsigned int get_ldx();
  /**
   * Returns the local size of the y dimension.
   */
  unsigned int get_ldy();
  /**
   * Returns the local size of the z dimension.
   */
  unsigned int get_ldz();


  /**
   * Returns the spacing in the x direction (mm)
   */
  delta_t get_deltax();
  /**
   * Returns the spacing in the y direction (mm)
   */
  delta_t get_deltay();
  /**
   * Returns the spacing in the z direction (mm)
   */
  delta_t get_deltaz();
  /**
   * Returns the spacing in time (s)
   */
  delta_t get_deltat();

  // Assign and get values of field components

  /**
   * Assign a value to Ex at some point in space.
   * @param x coordinate
   * @param y coordinate
   * @param z coordinate
   * @param value 
   */
  void set_ex(unsigned int x, unsigned int y, 
                     unsigned int z, field_t val);

  /**
   * Assign a value to Ey at some point in space.
   * @param x coordinate
   * @param y coordinate
   * @param z coordinate
   * @param value 
   */
  void set_ey(unsigned int x, unsigned int y, 
                     unsigned int z, field_t val);

  /**
   * Assign a value to Ez at some point in space.
   * @param x coordinate
   * @param y coordinate
   * @param z coordinate
   * @param value 
   */
  void set_ez(unsigned int x, unsigned int y, 
                     unsigned int z, field_t val);

  /**
   * Assign a value to Hx at some point in space.
   * @param x coordinate
   * @param y coordinate
   * @param z coordinate
   * @param value 
   */
  void set_hx(unsigned int x, unsigned int y, 
                     unsigned int z, field_t val);

  /**
   * Assign a value to Hy at some point in space.
   * @param x coordinate
   * @param y coordinate
   * @param z coordinate
   * @param value 
   */
  void set_hy(unsigned int x, unsigned int y, 
                     unsigned int z, field_t val);

  /**
   * Assign a value to Hz at some point in space.
   * @param x coordinate
   * @param y coordinate
   * @param z coordinate
   * @param value 
   */
  void set_hz(unsigned int x, unsigned int y, 
                     unsigned int z, field_t val);
};

#endif // GRID

