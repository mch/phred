
#include "SubdomainBc.hh"
#include "Grid.hh"

void SubdomainBc::apply(Face face, Grid &grid)
{
  region_t r = find_face(face, grid);
  MPI_Datatype t;

  t = grid.get_plane_dt(face);

  // Send away!
  send_recv(grid.get_face_start(face, EX, 1),
            grid.get_face_start(face, EX, 0), t);
  send_recv(grid.get_face_start(face, EY, 1),
            grid.get_face_start(face, EY, 0), t);
  send_recv(grid.get_face_start(face, EZ, 1),
            grid.get_face_start(face, EZ, 0), t);

  send_recv(grid.get_face_start(face, HX, 1),
            grid.get_face_start(face, HX, 0), t);
  send_recv(grid.get_face_start(face, HY, 1),
            grid.get_face_start(face, HY, 0), t);
  send_recv(grid.get_face_start(face, HZ, 1),
            grid.get_face_start(face, HZ, 0), t);
}

void SubdomainBc::send_recv(void *tx_ptr, void *rx_ptr, MPI_Datatype &t)
{
  MPI_Status status;

  if (rank_ > neighbour_) {
    MPI_Send(tx_ptr, 1, t, neighbour_, 1, MPI_COMM_WORLD);
    MPI_Recv(rx_ptr, 1, t, neighbour_, 1, MPI_COMM_WORLD, &status);
  } else {
    MPI_Recv(rx_ptr, 1, t, neighbour_, 1, MPI_COMM_WORLD, &status);
    MPI_Send(tx_ptr, 1, t, neighbour_, 1, MPI_COMM_WORLD);
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
