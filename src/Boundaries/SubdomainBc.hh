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

#ifndef SUBDOMAIN_BC_H
#define SUBDOMAIN_BC_H

#include "BoundaryCondition.hh"
#include "../Data.hh"

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
  void send_recv(const void *tx_ptr, void *rx_ptr, MPI_Datatype &t);

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
   * @param type the field to update
   */
  void apply(Face face, Grid &grid, FieldType type);

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

  /**
   * Returns a BoundayCondition type so the grid knows this is a subdomain.
   */
  virtual BoundaryCondition get_type() const;

};
#endif // SUBDOMAIN_BC_H
