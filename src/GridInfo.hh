#ifndef GRID_INFO_H
#define GRID_INFO_H

#include "Types.hh"

/**
 * This class is plain old data that is in class form "just in
 * case". Think of it as a structure. All members are public, but
 * this may change in the future. 
 *
 * This class holds data about grids. Number of cells in the total
 * grid, number of cells and starting point of the sub grid defined
 * for one processor, cell spacing in space, time step size, etc. 
 *
 * Information about the boundary conditions for the grids is also
 * held here. This allows the subdomaining algorithm to intellegently
 * reassign boundary conditions where required. 
 *
 * This object is more convenient for a parser to work with than the
 * full grid, and it makes the subdomain algorithm easier too. 
 */
class GridInfo
{
public:
  // Global grid size (i.e. all domains). These are *sizes*; the
  // maximum index into the {e,h}{x,y,z}_ arrays is these minus one. 
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

  // A grid is a cube with six faces. Those faces either need to have
  // boundary conditions, or they are subdomain boundaries and they
  // need to be shared with other processors. These arrays tell what
  // to do with each face. 
  //
  // 0 - Front (x = dimx, YZ plane)
  // 1 - Back (x = 0, YZ plane)
  // 2 - Left (y = 0, XZ plane)
  // 3 - Right (y = dimy, XZ plane)
  // 4 - Bottom (z = 0, XY plane)
  // 5 - Top (z = dimz, XY plane)
  //
  BoundaryCondition face_bc_[6]; // Boundary condition to apply
  int face_rank_[6]; // Rank of processor to talk to about this
		     // interface. 

  
  GridInfo() 
    : global_dimx_(0), global_dimy_(0), global_dimz_(0), 
      dimx_(0), dimy_(0), dimz_(0), 
      deltax_(0), deltay_(0), deltaz_(0), deltat_(0)
  {
    for (int i = 0; i < 6; i++) {
      face_bc_[i] = EWALL;
      face_rank_[i] = 0;
    }
  }

  ~GridInfo() {}

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
  inline void set_boundary(unsigned int face, BoundaryCondition bc)
  {
    face_bc_[face] = bc;
  }

  /**
   * Set the rank of the processor this face needs to be shared with. 
   *
   * @param face the face to assign the node rank to. One of FRONT, BACK,
   * LEFT, RIGHT, BOTTOM, TOP as defined in Types.hh
   *
   * @param rank the rank of the node to share this face with. 
   */ 
  inline void set_face_rank(unsigned int face, int rank)
  {
    face_rank_[face] = rank;
  }


};

#endif // GRID_INFO_H
