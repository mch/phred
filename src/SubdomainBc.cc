
#include "SubdomainBc.hh"
#include "Grid.hh"

void SubdomainBc::apply(Face face, Grid &grid)
{
  region_t r = find_face(face, grid);
  MPI_Datatype t;

  t = grid.get_plane_dt(face);

  // Send away!
  send_recv(grid.get_face_start(face, FC_EX, 1),
            grid.get_face_start(face, FC_EX, 0), t);
  send_recv(grid.get_face_start(face, FC_EY, 1),
            grid.get_face_start(face, FC_EY, 0), t);
  send_recv(grid.get_face_start(face, FC_EZ, 1),
            grid.get_face_start(face, FC_EZ, 0), t);

  send_recv(grid.get_face_start(face, FC_HX, 1),
            grid.get_face_start(face, FC_HX, 0), t);
  send_recv(grid.get_face_start(face, FC_HY, 1),
            grid.get_face_start(face, FC_HY, 0), t);
  send_recv(grid.get_face_start(face, FC_HZ, 1),
            grid.get_face_start(face, FC_HZ, 0), t);

  // Send and recieve and data we've be contracted to do. 
  vector<Data>::iterator tx_iter = tx_data_.begin();
  vector<Data>::iterator rx_iter = rx_data_.begin();
  vector<Data>::iterator tx_iter_end = tx_data_.end();

  while (tx_iter != tx_iter_end) 
  {
    send_recv((*tx_iter).get_ptr(), (*rx_iter).get_ptr(), 
              (*tx_iter).get_datatype());
  }
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

void SubdomainBc::add_tx_rx_pair(const Data &rx, const Data &tx)
{
  tx_data_.push_back(tx);
  rx_data_.push_back(rx);
}
