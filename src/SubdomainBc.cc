
#include "SubdomainBc.hh"
#include <iostream>

using namespace std;

void SubdomainBc::apply(Face face, Grid &grid)
{
  cout << "Supposed to apply subdomain boundary condition; send interface to "
       << neighbour_ << endl;

  region_t r = find_face(face, grid);
  MPI_Datatype t;
  void *start_ptr;

  switch (face)
  {
  case FRONT:
  case BACK:
    t = grid.get_yz_plane_dt();
    break;

  case LEFT:
  case RIGHT:
    t = grid.get_xz_plane_dt();
    break;

  case TOP:
  case BOTTOM:
    t = grid.get_xy_plane_dt();
    break;
  }

  // Send away!
  MPI_Status status;
  if (rank_ > neighbour_) {
    MPI_Send(start_ptr, 1, t, neighbour_, 1, MPI_COMM_WORLD);
    MPI_Recv(start_ptr, 1, t, neighbour_, 1, MPI_COMM_WORLD, &status);
  } else {
    // Temporary storage. Maybe a member var, so we don't have to alloc
    // and dealloc all the time? 
    int sz;
    MPI_Type_size(t, &sz);
    char *recv_buffer = new char[sz];
    
    MPI_Recv(recv_buffer, 1, t, neighbour_, 1, MPI_COMM_WORLD, &status);
    MPI_Send(start_ptr, 1, t, neighbour_, 1, MPI_COMM_WORLD);

    // Copy recieved data into grid

    // Cleanup
    delete[] recv_buffer;
  }

}

void SubdomainBc::set_neighbour(int neighbour)
{
  neighbour_ = neighbour;
}

int SubdomainBc::get_neighbour()
{
  return neighbour_;
}

void SubdomainBc::set_rank(int rank)
{
  rank_ = rank;
}

int SubdomainBc::get_rank()
{
  return rank_;
}
