#ifndef SUBDOMAIN_BC_H
#define SUBDOMAIN_BC_H

#include "BoundaryCondition.hh"

/**
 * This boundary condition talks to another rank and exchanges
 * information about the overlapping region with it. 
 * 
 * \bug IMPLEMENT ME!
 */
class SubdomainBc : public BoundaryCond
{
private:
protected:
  int neighbour_;

public:
  SubdomainBc() 
  {}

  ~SubdomainBc()
  {}

  /**
   * Exchanges the overlaping part of the grid with the adjacent processor. 
   *
   * @param face the face to apply the boundary condition to. 
   * @param grid the grid to apply the boundary condition to. 
   */
  void apply(Face face, Grid &grid);

  /**
   * Set the neighbour to talk to.
   *
   * @param neighbour the rank to talk to about this face.
   */
  void set_neighbour(int neighbour);

  /**
   * Returns the neighbour process rank for this face. 
   *
   * @return the neighbour rank
   */
  int get_neighbour();

};
#endif // SUBDOMAIN_BC_H
