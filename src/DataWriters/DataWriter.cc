/* 
   Phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004-2005 Matt Hughes <mhughe@uvic.ca>

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
{ 
  map<Variable *, Vardata *>::iterator iter = auxvardata_.begin();
  map<Variable *, Vardata *>::iterator iter_e = auxvardata_.end();

  for (; iter != iter_e; ++iter)
  {
    vector<MPI_Datatype>::iterator titer = (iter->second)->types_.begin();
    vector<MPI_Datatype>::iterator titer_e = (iter->second)->types_.end();
    
    for(; titer != titer_e; ++titer)
    {
      if (*titer != MPI_DATATYPE_NULL)
        MPI_Type_free(&(*titer));
    }

    if (iter->second->buffer_)
      delete (iter->second)->buffer_;

    delete (iter->second);
  }
}

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
  const Data &data = var.get_data();
  Vardata *auxvdata;

  MPI_Status status;

  map<Variable *, Vardata *>::iterator auxviter = auxvardata_.find(&var);

  if (auxvardata_.end() == auxviter)
  {
    auxvdata = new Vardata;
    auxvardata_[&var] = auxvdata;
    
    auxvdata->output_time_ = 0;
    auxvdata->types_ = gather_types(var);
    auxvdata->buffer_ = 0;
    auxvdata->buf_sz_ = 0;

    if (MPI_RANK == rank_)
    {
      int buffer_size = 1;
      const vector<Dimension> &dimensions = var.get_dimensions();
      for (vector<Dimension>::const_iterator iter = dimensions.begin();
           iter != dimensions.end(); ++iter)
        buffer_size = buffer_size * iter->global_len_;
      
      int type_sz = 0;
      MPI_Type_size(var.get_element_type(), &type_sz);

      auxvdata->buf_sz_ = buffer_size * type_sz;
      auxvdata->buffer_ = new char[auxvdata->buf_sz_];
    }
  } 
  else
  {
    auxvdata = auxviter->second;
  }

  // The number of items could potential change with every time step,
  // though usually it will be either 0 or some other constant number.
  vector<unsigned int> sizes = get_recieve_sizes(data);

  if (MPI_RANK != rank_)
  {
    if (data.get_num() > 0)
    {
      int sz;
      MPI_Type_size(data.get_datatype(), &sz);
      
      MPI_Send(const_cast<void *>(data.get_ptr()), 
               data.get_num(), data.get_datatype(), 
               rank_, 1, MPI_COMM_PHRED);
    } 
  }
  else 
  {
    unsigned int rcv_size = 0;

    for (vector<unsigned int>::iterator iter = sizes.begin();
         iter != sizes.end(); ++iter)
      rcv_size += *iter;

    if (rcv_size > 0) 
    {
      // Copy rank_ data into the buffer
      if (data.get_num() > 0) {

        if (MPI_SIZE == 1)
        {
          MPI_Sendrecv(const_cast<void *>(data.get_ptr()), 
                       data.get_num(), 
                       data.get_datatype(), 0, 1, 
                       static_cast<void *>(auxvdata->buffer_), 
                       sizes[0], MPI_CHAR, 0, 1,
                       MPI_COMM_PHRED, &status);
        } else {
          MPI_Sendrecv(const_cast<void *>(data.get_ptr()), 
                       data.get_num(), 
                       data.get_datatype(), 0, 1, 
                       static_cast<void *>(auxvdata->buffer_), 
                       1, auxvdata->types_[0], 0, 1,
                       MPI_COMM_PHRED, &status);
        }
      }
      
      // Recieve data from ranks =! rank_. This should really be
      // restricted to ranks in the current group etc...
      for (int i = 1; i < MPI_SIZE; i++)
      {
        if (sizes[i] > 0) {
          MPI_Recv(static_cast<void *>(auxvdata->buffer_), 
                   1, auxvdata->types_[i], 
                   i, 1, MPI_COMM_PHRED, &status);
        }
      }

      write_data(time_step, var, auxvdata->buffer_, 
                 auxvdata->buf_sz_);

      auxvdata->output_time_++;
    }
  }
}

vector<unsigned int> DataWriter::get_recieve_sizes(const Data &data)
{
  unsigned int nums_snd[2], *nums_recv;
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
             0, MPI_COMM_PHRED);

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

  for (; iter != iter_e; ++iter, ++dim_idx)
  {
    dim_lengths[dim_idx] = iter->local_len_;
    dim_starts[dim_idx] = iter->start_;
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
             num_dimensions, MPI_UNSIGNED, 0, MPI_COMM_PHRED);

  MPI_Gather(reinterpret_cast<void *>(dim_lengths), num_dimensions, 
             MPI_UNSIGNED, reinterpret_cast<void *>(recv_dim_lens), 
             num_dimensions, MPI_UNSIGNED, 0, MPI_COMM_PHRED);


  if (MPI_RANK == rank_)
  {
    // Construct data types for each node's incoming data. 
    int *dim_sizes, *dim_sub_sizes, *starts;
    int *displacements = 0;
    unsigned int num_displacements = 0;

    dim_sizes = new int[num_dimensions];
    dim_sub_sizes = new int[num_dimensions];
    starts = new int[num_dimensions];

    for (dim_idx = 0; dim_idx < num_dimensions; dim_idx++)
    {
      dim_sizes[dim_idx] = dimensions[dim_idx].global_len_;
    }

    for (unsigned int node_idx = 0; node_idx < MPI_SIZE; node_idx++)
    {
      MPI_Datatype arr_type;

      int sz = 1;
      for (unsigned int node_dim_idx = 0, dim_idx = node_idx * num_dimensions; 
           dim_idx < (node_idx + 1) * num_dimensions; 
           dim_idx++, node_dim_idx++)
      {
        dim_sub_sizes[node_dim_idx] = recv_dim_lens[dim_idx];
        starts[node_dim_idx] = recv_dim_starts[dim_idx];

        sz = sz * dim_sub_sizes[node_dim_idx];
      }

      if (sz > 0)
      {
        MPI_Type_create_subarray(num_dimensions, dim_sizes, 
                                 dim_sub_sizes, starts, 
                                 MPI_ORDER_C, var.get_element_type(), 
                                 &arr_type);
        MPI_Type_commit(&arr_type);    
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
