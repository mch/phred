#ifndef SUBDOMAIN_BC_H
#define SUBDOMAIN_BC_H

#include "BoundaryCondition.hh"
#include "Data.hh"

#include <vector>

using namespace std;

/**
 * This boundary condition talks to another rank and exchanges
 * information about the overlapping region with it. 
 * 
 */
class SubdomainBc : public BoundaryCond
{
private:
protected:
  int neighbour_;
  int rank_;

  /**
   * A list of items to transmit each time the boundary condition is
   * applied. The communication is assumed to be bidirectional and
   * non-overlapping.
   */ 
  vector<RxTxData> rx_tx_data_;

  /**
   * A helper function that sends arrays of data between two ranks. 
   * @param tx_ptr pointer to data to send
   * @param tx_ptr pointer to area to recieve to 
   * @param dt MPI Derived data type describing the data
   */
  void send_recv(void *tx_ptr, void *rx_ptr, MPI_Datatype &t);

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

  /**
   * Set our own rank.
   *
   * @param rank our own rank
   */
  void set_rank(int rank);

  /**
   * Returns out rank. 
   *
   * @return our rank
   */
  int get_rank();

  /**
   * Add a tx/rx data pair to be exchanged each time this boundary
   * condition is applied. The datatype used by both Data objects MUST
   * be the same.  
   *
   * @param x A RxTxData object describing the data we will recieve.
   */
  void add_tx_rx_data(const RxTxData &x);

};
#endif // SUBDOMAIN_BC_H
