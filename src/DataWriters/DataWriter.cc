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
  unsigned int nums_snd[2], *nums_recv;
  char *ptr, *ptr_head; 
  int sz; 
  const Data &data = var.get_data();

  MPI_Status status;

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

  if (rank_ != 0)
  {
    if (data.get_num() > 0)
      MPI_Send(data.get_ptr(), data.get_num(), data.get_datatype(), 
               0, 1, MPI_COMM_WORLD);
  } 
  else 
  {
    unsigned int total = 0;

    // Total number of bytes to recieve...
    for (int i = 0; i < size_; i++)
      total += nums_recv[i * 2] * nums_recv[(i*2)+1];

    if (total > 0) 
    {
    
      // A buffer to hold the result
      ptr_head = ptr = new char[total];
      
      // Copy rank 0 data into the buffer
      if (data.get_num() > 0) {
        MPI_Sendrecv(data.get_ptr(), data.get_num(), 
                     data.get_datatype(), 0, 1, 
                     static_cast<void *>(ptr), 
                     nums_recv[0] * nums_recv[1], MPI_CHAR, 0, 1,
                     MPI_COMM_WORLD, &status);
        ptr += nums_recv[0] * nums_recv[1];
      }
      
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

  if (rank_ == 0)
    delete[] nums_recv;  
}
