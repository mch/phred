#ifndef GRID_INFO_H
#define GRID_INFO_H

#include "Types.hh"
#include "BoundaryCondition.hh"
#include "Ewall.hh"
#include "Hwall.hh"
#include "SubdomainBc.hh"
#include "Pml.hh"

#include "counted_ptr.hh"

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
protected:
  // THIS MAY BE A PRIME SMART POINTER CANDIDATE, due to the copying complexity. 
  counted_ptr<BoundaryCond> *face_bc_[6]; // Boundary condition to apply

  /**
   * GridInfo will occasionally have to take responsibility for
   * deleting boundary conditions that get left with us, such as those
   * created by the domain decomposition algorithm. This array tells
   * us whats what (true is GridInfo owns the bc).
   */
  //bool face_ptr_owner_[6];
  
  
public:
  GridInfo();

  /**
   * Copy constructor, to properly handle the dynamically allocated
   * boundary conditions. 
   *
   * @param info the GridInfo object to be copied. 
   */
  GridInfo(const GridInfo &info);

  ~GridInfo();

  /**
   * Assignment operator. To handle the dynamically allocated
   * boundary conditions properly.
   */
  GridInfo& operator=(const GridInfo &info);

  /**
   * Set the boundary condition on one of the faces of this grid. 
   *
   * @param face the face to assign the boundary to. One of FRONT, BACK,
   * LEFT, RIGHT, BOTTOM, TOP as defined in Types.hh
   *
   * @param bc the boundary condition to apply. 
   *
   * @return a BoundaryCond object of the type required, in which the
   * specifics of the boundary condition can be stored.
   */ 
  void set_boundary(Face face, BoundaryCond *bc, bool take_ownership = false);

  /**
   * Set a PML boundary
   * @param face the face to apply the pml to 
   * @param thickness the number of cells thick the pml is
   * @param var the pml variation profile
   * @param nrml_refl amount of normal reflection (try 1.0)
   * @return a point to the pml object if you want to mess with it some more. 
   */
  Pml *set_pml_boundary(Face face, unsigned int thickness, 
                        PmlVariation_t var, float nrml_refl);

  /**
   * Returns the type of boundary assigned to a face.
   *
   * @return BoundaryCondition from Types.hh
   */
  inline const BoundaryCondition get_bc_type(Face face)
  {
    return face_bc_[face]->get()->get_type();
  }

  /** 
   * Returns a reference to the boundary condition object for a face. 
   *
   * @return ref to a BoundaryCond
   */
  inline BoundaryCond& get_boundary(Face face)
  {
    return *(face_bc_[face]->get());
  }

  /**
   * Returns the face thickness for a boundary condition
   *
   * @return an unsigned int, the thickness of the boundary condition.
   */
  unsigned int get_face_thickness(Face face);

  /**
   * Apply the boundary conditions to the grid. 
   *
   * @param the grid to apply to 
   */
  void apply_boundaries(Grid &grid, FieldType type);

};

#endif // GRID_INFO_H
