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

#include "DataWriter.hh"
#include "../Globals.hh"

void DataWriter::set_filename(string filename)
{
  filename_ = filename;
}

void DataWriter::handle_data(unsigned int time_step, 
                             map<string, Variable *> &data)
{
  map<string, Variable *>::iterator iter = data.begin();
  map<string, Variable *>::iterator iter_e = data.end();

  for (iter = data.begin(); iter != iter_e; ++iter)
    gather_data(time_step, *(*iter).second);
}

void DataWriter::gather_data(unsigned int time_step, Variable &var)
{
  // Is the data on the right rank? 
  char *ptr, *ptr_head; 
  int buffer_size, sz; 
  const Data &data = var.get_data();

  MPI_Status status;


  if (rank_ != 0)
  {
    if (data.get_num() > 0)
      MPI_Send(data.get_ptr(), data.get_num(), data.get_datatype(), 
               0, 1, MPI_COMM_WORLD);
  } 

  else 
  {
    unsigned int rcv_size = 0;
    vector<unsigned int> rcv_bytes = get_recieve_sizes(data);


    if (rcv_size > 0) 
    {
      vector<pair<unsigned int, unsigned int> > sizes = gather_sizes(var);

      // A buffer to hold the result
      ptr_head = ptr = new char[buffer_size];
      
      // Copy rank 0 data into the buffer
      // THIS IS BUGGY FOR MPI_SIZE > 1!!!! It does not take size and
      // start of region into account. 
      if (data.get_num() > 0) {
        MPI_Sendrecv(data.get_ptr(), data.get_num(), 
                     data.get_datatype(), 0, 1, 
                     static_cast<void *>(ptr), 
                     nums_recv[0] * nums_recv[1], MPI_CHAR, 0, 1,
                     MPI_COMM_WORLD, &status);
        ptr += nums_recv[0] * nums_recv[1];
      }
      
      // Recieve data from ranks > 0 
      // THIS IS ALSO BUGGY
      for (int i = 1; i < size_; i++)
      {
        if (nums_recv[i*2] > 0) {
          MPI_Recv(ptr, nums_recv[i*2] * nums_recv[i*2+1], MPI_CHAR, 
                   i, 1, MPI_COMM_WORLD, &status);
          ptr += nums_recv[i*2] * nums_recv[i*2+1];
        }
      }
      
      MPI_Datatype t = data.get_datatype();
        
      write_data(time_step, var, t, ptr_head, total);

      delete[] ptr_head;
    }
  }

  delete[] dim_starts;
  delete[] dim_lengths;

  if (rank_ == 0)
  {
    delete[] recv_dim_starts;
    delete[] recv_dim_lengths;
  }
}

vector<unsigned int> DataWriter::get_recieve_sizes(const Data &data)
{
  unsigned int nums_snd[2], *nums_recv;
  unsigned int total = 0;
  vector<unsigned int> ret;

  if (rank_ == 0)
    nums_recv = new unsigned int[size_ * 2];

  // Every process sends the number of items it has to contribute to
  // the result to rank 0. 
  nums_snd[0] = data.get_num();
  MPI_Type_size(data.get_datatype(), &sz);
  nums_snd[1] = static_cast<unsigned int>(sz);
  
  MPI_Gather(static_cast<void *>(&nums_snd), 2, MPI_UNSIGNED, 
             static_cast<void *>(nums_recv), 2, MPI_UNSIGNED, 
             0, MPI_COMM_WORLD);

  if (rank_ == 0)
  {
    for (int i = 0; i < size_; i++)
      ret.push_back(nums_recv[i * 2] * nums_recv[(i*2)+1]);
    
    delete[] nums_recv;
  }

  return ret;
}

vector<pair<unsigned int, unsigned int> > gather_sizes(const Variable &var)
{
  // Every process must let us know where the data it is sending
  // starts, and what it's dimensions are, so we can fit it into the
  // right place in the buffer. These numbers are in terms of the 
  const vector<Dimension> dimensions = var.get_dimensions();
  vector<Dimension>::const_iterator iter = dimensions.begin();
  vector<Dimension>::const_iterator iter_e = dimensions.end();
  unsigned int num_dimensions = dimensions.size();
  unsigned int *dim_starts = new unsigned int[num_dimensions];
  unsigned int *dim_lengths = new unsigned int[num_dimensions];
  unsigned int dim_idx = 0;
  vector<pair<unsigned int, unsigned int> > ret;

  cerr << "Rank " << MPI_RANK << " is writing a variable with "
       << num_dimensions << " dimensions. \n\tstart\tlength:\n";

  for (; iter != iter_e; ++iter, ++dim_idx)
  {
    if (data.get_num() > 0)
    {
      dim_lengths[dim_idx] = iter->local_len_;
      dim_starts[dim_idx] = iter->start_;
    } else {
      dim_lengths[dim_idx] = 0;
      dim_starts[dim_idx] = 0;
    }
    cerr << "\t" << dim_starts[dim_idx] << "\t" << dim_lengths[dim_idx] << "\n";
  }

  if (rank_ == 0)
  {
    unsigned int *recv_dim_starts = new unsigned int[num_dimensions * MPI_SIZE];
    unsigned int *recv_dim_lens = new unsigned int[num_dimensions * MPI_SIZE];
  }

  MPI_Gather(reinterpret_cast<void *>(dim_starts), num_dimensions, 
             MPI_UNSIGNED, reinterpret_cast<void *>(recv_dim_starts), 
             num_dimensions, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  MPI_Gather(reinterpret_cast<void *>(dim_lengths), num_dimensions, 
             MPI_UNSIGNED, reinterpret_cast<void *>(recv_dim_lens), 
             num_dimensions, MPI_UNSIGNED, 0, MPI_COMM_WORLD);


  if (MPI_RANK == 0)
  {
    unsigned int *total_dim_lens = new unsigned int[num_dimensions];
    memset(total_dim_lens, 0, num_dimensions * sizeof(unsigned int));
    
    cerr << "Rank 0 has gathered all data to it:\n\tstart\tlength:\n";
    for (dim_idx = 0; dim_idx < num_dimensions * MPI_SIZE; dim_idx++)
    {
      cerr << "\t" << recv_dim_starts[dim_idx] << "\t" 
           << recv_dim_lens[dim_idx] << "\n";
      unsigned int offset = dim_idx % num_dimensions;
      total_dim_lens[offset] += recv_dim_lens[dim_idx];
    }

    cerr << "Totalled dimension lengths... \n";
    for (dim_idx = 0; dim_idx < num_dimensions; dim_idx++)
      cerr << "\t" << total_dim_lens[dim_idx] << "\n";
  }

  // DO ANY REGIONS OVERLAP? Is the total number of items to be
  // recieved, as far as the lengths is concerned, deviate from
  // expected?

}
