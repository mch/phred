/* 
   Phred - Phred is a parallel finite difference time domain
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
#include "../Contiguous.hh"

DataWriter::DataWriter() 
  : rank_(0)
{ }

DataWriter::~DataWriter()
{ }

void DataWriter::add_scalar(const char *name, double &value)
{
  throw DataWriterException("This data writer does not support add_scalar.");
}

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

  // How much data will we recieve? 
  vector<unsigned int> rcv_bytes = get_recieve_sizes(data);

  if (!var.have_node_dtypes_)
  {
    var.node_types_ = gather_types(var);
    var.have_node_dtypes_ = true;
  }

//   cerr << "Rank " << MPI_RANK << " has " << data.get_num() 
//        << " items to contribute to variable " << var.get_name() << endl;

  if (MPI_RANK != rank_)
  {
    if (data.get_num() > 0)
    {
      int sz;
      MPI_Type_size(data.get_datatype(), &sz);
      
//       cerr << "Rank " << MPI_RANK << " is sending " << sz << " bytes..."
//            << endl;

      MPI_Send(data.get_ptr(), data.get_num(), data.get_datatype(), 
               0, 1, MPI_COMM_WORLD);

//       cerr << "Rank " << MPI_RANK << " sent " << sz << " bytes..."
//            << endl;      
    } 
//     else
//     {
//       cerr << "Rank " << MPI_RANK << " has no data for variable " 
//            << var.get_name() << endl;
//     }
  } 
  else 
  {
    unsigned int rcv_size = 0;

    for (vector<unsigned int>::iterator iter = rcv_bytes.begin();
         iter != rcv_bytes.end(); ++iter)
      rcv_size += *iter;

#ifdef DEBUG
//     cerr << "Datawriter on rank 0 will recieve a total of " 
//          << rcv_size << " bytes. \n";
#endif

    if (rcv_size > 0) 
    {
      // Construct the datatype that the subclasses will use to
      // interpert the buffered data. 
      MPI_Datatype global_type = data.get_global_datatype();

      if (MPI_DATATYPE_NULL == global_type)
      {
        const vector<Dimension> &dimensions = var.get_dimensions();

        unsigned int buffer_size = 1;
        for (vector<Dimension>::const_iterator iter = dimensions.begin();
             iter != dimensions.end(); ++iter)
          buffer_size = buffer_size * iter->global_len_;

//         cerr << "CONTIGUOUS BUFFER IS " << buffer_size
//              << " elements long!" << endl;

        MPI_Type_contiguous(buffer_size, var.get_element_type(),
                            &global_type);
        MPI_Type_commit(&global_type);
      }
      
      int buffer_size = 0;
      MPI_Type_size(global_type, &buffer_size);

      if (buffer_size < rcv_size)
        throw DataWriterException("Global data type size is smaller\
 than the amount of data to be recieved!");
      
      ptr_head = ptr = new char[buffer_size];
      
//       cerr << "Rank " << rank_ << " reciever buffer size is " 
//            << buffer_size << endl;

      // Copy rank_ data into the buffer
      if (data.get_num() > 0) {

        int sz = 0;
        MPI_Type_size(data.get_datatype(), &sz);
//         cerr << "Rank " << rank_ << " is swizilling " 
//              << data.get_num() << " items (" << sz 
//              << " bytes) with itself. Buffer is " 
//              << rcv_bytes[0] << " bytes. " << endl;

        if (MPI_SIZE == 1)
        {
          MPI_Sendrecv(data.get_ptr(), data.get_num(), 
                       data.get_datatype(), 0, 1, 
                       static_cast<void *>(ptr), 
                       rcv_bytes[0], MPI_CHAR, 0, 1,
                       MPI_COMM_WORLD, &status);
        } else {
          int sz = 0;
          MPI_Type_size(var.node_types_[0], &sz);

//           cerr << "MPI type is " << sz << " bytes..." << endl;

          MPI_Sendrecv(data.get_ptr(), data.get_num(), 
                       data.get_datatype(), 0, 1, 
                       static_cast<void *>(ptr), 
                       1, var.node_types_[0], 0, 1,
                       MPI_COMM_WORLD, &status);

//           cerr << "Buffered " << rcv_bytes[0] << " from rank 0. " 
//                << "Or is it " << sz << " bytes?" << endl;
        }
      }
      
      // Recieve data from ranks =! rank_
      for (int i = 1; i < MPI_SIZE; i++)
      {
        if (rcv_bytes[i] > 0) {
          // Construct a derived data type that places the recieved
          // data in the correct location inside the buffer. 
          int sz = 0;
          MPI_Type_size(var.node_types_[i], &sz);
//           cerr << "Rank 0 is recieving " << rcv_bytes[i] 
//                << " (" << sz << "?) bytes from rank " << i << endl;

          MPI_Recv(static_cast<void *>(ptr), 1, var.node_types_[i], 
                   i, 1, MPI_COMM_WORLD, &status);

//           cerr << "Rank 0 recieved " << rcv_bytes[i] 
//                << " (" << sz << "?) bytes from rank " << i << endl;

        }
      }

//       cout << "TEST -------------------\n";
//       for (int i=0; i < 2*MPI_SIZE; i++)
//       {
//         for (int j = 0; j < 8; j++)
//           cout << "\t" << reinterpret_cast<float *>(ptr)[i + j*2];
//         cout << endl;
//       }
//       cout << "-------------------\n";
      
      write_data(time_step, var, global_type, ptr_head, buffer_size);
      var.output_time_++;

      MPI_Type_free(&global_type);

      delete[] ptr_head;
    }
  }

}

vector<unsigned int> DataWriter::get_recieve_sizes(const Data &data)
{
  unsigned int nums_snd[2], *nums_recv;
  unsigned int total = 0;
  vector<unsigned int> ret;
  int sz;

  if (MPI_RANK == rank_)
    nums_recv = new unsigned int[MPI_SIZE * 2];

  // Every process sends the number of items it has to contribute to
  // the result to rank 0. 
  nums_snd[0] = data.get_num();
  MPI_Type_size(data.get_datatype(), &sz);
  nums_snd[1] = static_cast<unsigned int>(sz);
  
  MPI_Gather(static_cast<void *>(&nums_snd), 2, MPI_UNSIGNED, 
             static_cast<void *>(nums_recv), 2, MPI_UNSIGNED, 
             0, MPI_COMM_WORLD);

  if (MPI_RANK == rank_)
  {
    for (int i = 0; i < MPI_SIZE; i++)
      ret.push_back(nums_recv[i * 2] * nums_recv[(i*2)+1]);
    
    delete[] nums_recv;
  }

  return ret;
}

vector<MPI_Datatype> DataWriter::gather_types(const Variable &var)
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
  vector<MPI_Datatype> ret;
  const Data &data = var.get_data();

//   cerr << "Rank " << MPI_RANK << " is writing the variable "
//        << var.get_name() << " with "
//        << num_dimensions << " dimensions. \n\tstart\tlength:\n";

  for (; iter != iter_e; ++iter, ++dim_idx)
  {
    dim_lengths[dim_idx] = iter->local_len_;
    dim_starts[dim_idx] = iter->start_;

//     cerr << "\t" << dim_starts[dim_idx] << "\t" 
//          << dim_lengths[dim_idx] << "\n";
  }

  unsigned int *recv_dim_starts = 0;
  unsigned int *recv_dim_lens = 0;

  if (MPI_RANK == rank_)
  {
    recv_dim_starts = new unsigned int[num_dimensions * MPI_SIZE];
    recv_dim_lens = new unsigned int[num_dimensions * MPI_SIZE];
  }

  MPI_Gather(reinterpret_cast<void *>(dim_starts), num_dimensions, 
             MPI_UNSIGNED, reinterpret_cast<void *>(recv_dim_starts), 
             num_dimensions, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  MPI_Gather(reinterpret_cast<void *>(dim_lengths), num_dimensions, 
             MPI_UNSIGNED, reinterpret_cast<void *>(recv_dim_lens), 
             num_dimensions, MPI_UNSIGNED, 0, MPI_COMM_WORLD);


  if (MPI_RANK == rank_)
  {
    // Construct data types for each node's incoming data. 
    unsigned int *dim_sizes, *dim_sub_sizes, *starts;
    int *displacements = 0;
    unsigned int num_displacements = 0;

    dim_sizes = new unsigned int[num_dimensions];
    dim_sub_sizes = new unsigned int[num_dimensions];
    starts = new unsigned int[num_dimensions];

//     cerr << "Total N-dimensional array size:" << endl;
    for (dim_idx = 0; dim_idx < num_dimensions; dim_idx++)
    {
      dim_sizes[dim_idx] = dimensions[dim_idx].global_len_;
//       cerr << "\t" << dim_sizes[dim_idx];
    }
//     cerr << endl << endl;

    for (unsigned int node_idx = 0; node_idx < MPI_SIZE; node_idx++)
    {
      MPI_Datatype arr_type;

//       cerr << "Setting up data type for data from node " << node_idx
//            << endl;

      int sz = 1;
      for (unsigned int node_dim_idx = 0, dim_idx = node_idx * num_dimensions; 
           dim_idx < (node_idx + 1) * num_dimensions; 
           dim_idx++, node_dim_idx++)
      {
        dim_sub_sizes[node_dim_idx] = recv_dim_lens[dim_idx];
        starts[node_dim_idx] = recv_dim_starts[dim_idx];

        sz = sz * dim_sub_sizes[node_dim_idx];

//         cerr << "Dimension " << node_dim_idx << " starts at "
//              << starts[node_dim_idx] << " and is " 
//              << dim_sub_sizes[node_dim_idx] << " long." << endl;
      }

      Contiguous::compute_displacements(num_dimensions, dim_sizes, starts, 
                                        dim_sub_sizes, &num_displacements, 
                                        &displacements);

//       cout << "Computed displacements, got " 
//            << num_displacements << " back, each has a length of "
//            << dim_sub_sizes[num_dimensions - 1] << ":\n";
      
//       for (int i = 0; i < num_displacements; i++)
//         cout << "\t" << displacements[i];

//       cout << endl << endl;
      
      if (sz > 0)
      {
        int *lengths = new int[num_displacements];

        for (int i = 0; i < num_displacements; i++)
          lengths[i] = dim_sub_sizes[num_dimensions - 1];

        MPI_Type_indexed(num_displacements, lengths, displacements, 
                         var.get_element_type(), &arr_type);

        MPI_Type_commit(&arr_type);

        delete[] lengths;
      }
      else
        arr_type = MPI_DATATYPE_NULL;

      ret.push_back(arr_type);

      delete[] displacements;

    }
    
    delete[] dim_sizes;
    delete[] dim_sub_sizes;
    delete[] starts;

    delete[] recv_dim_starts;
    delete[] recv_dim_lens;
  }

  // DO ANY REGIONS OVERLAP? Is the total number of items to be
  // recieved, as far as the lengths is concerned, deviate from
  // expected? Does it matter? 

  delete[] dim_starts;
  delete[] dim_lengths;

  return ret;
}

ostream& DataWriter::to_string(ostream &os) const
{
  return os << "A DataWriter of indeterminate type, writing to '"
            << filename_ << "' on rank " << rank_;
}
