#include "DataWriter.hh"

void DataWriter::handle_data(unsigned int time_step, Data &data)
{
  gather_data(data);
}

void DataWriter::gather_data(Data &data)
{
  // Is the data on the right rank? 
  unsigned int nums_snd[2], *nums_recv;
  char *ptr, *ptr_head; 
  int sz; 

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
        
      write_data(data, t, ptr_head, total);

      delete[] ptr_head;
    }
  }

  if (rank_ == 0)
    delete[] nums_recv;  
}
