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

  if (num_dims > 1)
    throw std::exception(); // "No more than one dimension allowed!");
  
  if (num_dims == 0)
    throw std::exception(); //"Result must have exactly one dimension.");
  
  dim_len_ = result.get_dim_lengths()[0];

}

void AsciiDataWriter::handle_data(unsigned int time_step, Data &data)
{
  // Is the data on the right rank? 
  unsigned int nums_snd[2], *nums_recv;
  void *ptr; 
  
  MPI_Status status;

  if (rank_ == 0)
    nums_recv = new unsigned int[size_ * 2];

  // Every process sends the number of items it has to contribute to
  // the result to rank 0. 
  nums_snd[0] = data.get_num();
  nums_snd[1] = data.get_num_bytes();

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

    for (int i = 0; i < size_; i++)
      total += nums_recv[i * 2] * nums_recv[(i*2)+1];

    ptr = new char[total];

    memcpy(ptr, data.get_ptr(), data.get_num() * data.get_num_bytes());
    (static_cast<char *>(ptr)) += data.get_num() * data.get_num_bytes();

    for (int i = 1; i < size_; i++)
    {
      if (nums_recv[i*2] > 0)
        MPI_Recv(ptr, nums_recv[i*2], data.get_datatype(), 
                 i, 1, MPI_COMM_WORLD, &status);
    }

    if (!file_.is_open()) 
    {
      init();
      
      if (!file_.is_open()) 
        throw std::exception(); // "Unable to open output file!");
    }

    MPI_Datatype t = data.get_datatype();

    write_data(t, ptr);
    file_ << endl;

    delete[] static_cast<char *>(ptr);
  }

  if (rank_ == 0)
    delete[] nums_recv;
}

void *AsciiDataWriter::write_data(MPI_Datatype t, void *ptr)
{
  if (t == MPI_CHAR) 
    file_ << *(static_cast<signed char *>(ptr)++) << "\t";
  else if (t == MPI_SHORT)
    file_ << *(static_cast<signed short int *>(ptr)++) << "\t";
  else if (t == MPI_INT)
    file_ << *(static_cast<signed int *>(ptr)++) << "\t";
  else if (t == MPI_LONG)
    file_ << *(static_cast<signed long int *>(ptr)++) << "\t";
  else if (t == MPI_UNSIGNED_CHAR)
    file_ << *(static_cast<unsigned char *>(ptr)++) << "\t";
  else if (t == MPI_UNSIGNED_SHORT)
    file_ << *(static_cast<unsigned short int *>(ptr)++) << "\t";
  else if (t == MPI_UNSIGNED)
    file_ << *(static_cast<unsigned int *>(ptr)++) << "\t";
  else if (t == MPI_UNSIGNED_LONG)
    file_ << *(static_cast<unsigned long int *>(ptr)++) << "\t";
  else if (t == MPI_FLOAT)
    file_ << *(static_cast<float *>(ptr)++) << "\t";
  else if (t == MPI_DOUBLE)
    file_ << *(static_cast<double *>(ptr)++) << "\t";
  else if (t == MPI_LONG_DOUBLE)
    file_ << *(static_cast<long double *>(ptr)++) << "\t";
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
      for (int j = 0; j < ints[i]; j++)
        ptr = write_data(dts[i], ptr);
      i++;
    }

    delete[] ints;
    delete[] aints;
    delete[] dts;
  }
  
  return ptr;
}
