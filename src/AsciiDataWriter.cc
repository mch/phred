#include "AsciiDataWriter.hh"
#include <string.h>

AsciiDataWriter::AsciiDataWriter(int rank, int size)
  : DataWriter(rank, size)
{}

AsciiDataWriter::~AsciiDataWriter()
{
  if (file_.is_open())
    file_.close();
}

void AsciiDataWriter::init()
{
  if (filename_.length() > 0 && rank_ == 0)
  {
    file_.open(filename_.c_str(), ofstream::out);
    file_.flags(ios_base::scientific);
  }
}

void AsciiDataWriter::deinit()
{

}

void AsciiDataWriter::add_variable(Result &result)
{
  unsigned int num_dims = result.get_num_dims();

  if (num_dims == 0)
    throw std::exception(); //"Result must have at least one dimension.");
  
  dim_len_ = 0;
  for (int i = 0; i < num_dims; i++)
    dim_len_ += result.get_dim_lengths()[i];

}

void AsciiDataWriter::handle_data(unsigned int time_step, Data &data)
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
                     static_cast<void *>(ptr), nums_recv[1], MPI_CHAR, 0, 1,
                     MPI_COMM_WORLD, &status);
        ptr += sz;
      }
      
      for (int i = 1; i < size_; i++)
      {
        if (nums_recv[i*2] > 0) {
          MPI_Recv(ptr, nums_recv[(i+1)*2], MPI_CHAR, 
                   i, 1, MPI_COMM_WORLD, &status);
          ptr += nums_recv[(i+1)*2];
        }
      }
      
      if (!file_.is_open()) 
      {
        init();
        
        if (!file_.is_open()) 
          throw std::exception(); // "Unable to open output file!");
      }
      
      MPI_Datatype t = data.get_datatype();
      
      write_data(t, ptr_head, total);
      file_ << endl;
      
      delete[] ptr_head;
    }
  }

  if (rank_ == 0)
    delete[] nums_recv;
}

void *AsciiDataWriter::write_data(MPI_Datatype t, void *ptr, 
                                  unsigned int len)
{
  if (t == MPI_CHAR) 
    ptr = write_array<signed char>(ptr, len);
  else if (t == MPI_SHORT)
    ptr = write_array<signed short int>(ptr, len);
  else if (t == MPI_INT)
    ptr = write_array<signed int>(ptr, len);
  else if (t == MPI_LONG)
    ptr = write_array<signed long int>(ptr, len);
  else if (t == MPI_UNSIGNED_CHAR)
    ptr = write_array<unsigned char>(ptr, len);
  else if (t == MPI_UNSIGNED_SHORT)
    ptr = write_array<unsigned short int>(ptr, len);
  else if (t == MPI_UNSIGNED)
    ptr = write_array<unsigned int>(ptr, len);
  else if (t == MPI_UNSIGNED_LONG)
    ptr = write_array<unsigned long int>(ptr, len);
  else if (t == MPI_FLOAT)
    ptr = write_array<float>(ptr, len);
  else if (t == MPI_DOUBLE)
    ptr = write_array<double>(ptr, len);
  else if (t == MPI_LONG_DOUBLE)
    ptr = write_array<long double>(ptr, len);
  else if (t == MPI_PACKED)
    file_ << "packed\t";
  else 
  {
    int num_ints, num_addrs, num_dts, combiner; 

    MPI_Type_get_envelope(t, &num_ints, &num_addrs, &num_dts, 
                          &combiner);

    int *ints;
    MPI_Aint *aints;
    MPI_Datatype *dts;

    ints = new int[num_ints];
    aints = new MPI_Aint[num_addrs];
    dts = new MPI_Datatype[num_dts];

    MPI_Type_get_contents(t, num_ints, num_addrs, num_dts, 
                          ints, aints, dts);

    int i = 0;
    while (i < num_dts && dts[i] > 0)
    {
      ptr = write_data(dts[i], ptr, ints[i]);
      i++;
    }

    delete[] ints;
    delete[] aints;
    delete[] dts;
  }
  
  return ptr;
}

template <class T>
void *AsciiDataWriter::write_array(void *ptr, unsigned int len)
{
  T *r = static_cast<T *>(ptr);

  for (unsigned int i = 0; i < len; i++)
  {
    file_ << *(r++) << "\t";
  }

  return static_cast<void *>(r);
}
