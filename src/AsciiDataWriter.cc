#include "AsciiDataWriter.hh"
#include <string.h>
#include <mpi.h>

AsciiDataWriter::AsciiDataWriter(int rank, int size)
  : DataWriter(rank, size)
{}

AsciiDataWriter::~AsciiDataWriter()
{
  deinit();
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
  if (rank_ == 0 && file_.is_open())
    file_.close();
}

void AsciiDataWriter::add_variable(Result &result)
{
  const vector<int> &dim_lens = result.get_dim_lengths();

  if (dim_lens.size() == 0)
    throw std::exception(); //"Result must have at least one dimension.");
  
  dim_len_ = 0;
  for (int i = 0; i < dim_lens.size(); i++)
    dim_len_ += dim_lens[i];

}

void *AsciiDataWriter::write_data(Data &data, MPI_Datatype t, 
                                  void *ptr, unsigned int len)
{
  if (!file_.is_open()) 
  {
    init();
    
    if (!file_.is_open()) 
      throw std::exception(); // "Unable to open output file!");
  }

  ptr = write_data(t, ptr, len);
  file_ << endl;

  return ptr;
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
