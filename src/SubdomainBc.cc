
#include "SubdomainBc.hh"
#include "Grid.hh"
#include <iostream>

using namespace std;

void SubdomainBc::apply(Face face, Grid &grid)
{
  cout << "Supposed to apply subdomain boundary condition; send interface to "
       << neighbour_ << endl;

  region_t r = find_face(face, grid);
  MPI_Datatype t;

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
  send_recv(grid.get_face_start(face, EX), t);
  send_recv(grid.get_face_start(face, EY), t);
  send_recv(grid.get_face_start(face, EZ), t);
  send_recv(grid.get_face_start(face, HX), t);
  send_recv(grid.get_face_start(face, HY), t);
  send_recv(grid.get_face_start(face, HZ), t);
}

void SubdomainBc::send_recv(void *ptr, MPI_Datatype &t)
{
  MPI_Status status;

  if (rank_ > neighbour_) {
    MPI_Send(ptr, 1, t, neighbour_, 1, MPI_COMM_WORLD);
    MPI_Recv(ptr, 1, t, neighbour_, 1, MPI_COMM_WORLD, &status);
  } else {
    // Temporary storage. Maybe a member var, so we don't have to alloc
    // and dealloc all the time? 
    int sz, num_items;
    MPI_Type_size(t, &sz);

    num_items = sz / sizeof(field_t);

    // Create a contigous type to recieve buffered data with. 
    MPI_Datatype recv_t; 
    
    MPI_Type_contiguous(num_items, GRID_MPI_TYPE, &recv_t);
    MPI_Type_commit(&recv_t);

    field_t *recv_buffer = new field_t[num_items];
    
    MPI_Recv(recv_buffer, 1, recv_t, neighbour_, 1, MPI_COMM_WORLD, &status);
    MPI_Send(ptr, 1, t, neighbour_, 1, MPI_COMM_WORLD);

    // Copy recieved data into grid
    int pos = 0;
    MPI_Unpack(recv_buffer, num_items, &pos, ptr, 1, t, MPI_COMM_WORLD);

    // Cleanup
    delete[] recv_buffer;

    MPI_Type_free(&recv_t); // !?
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
