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

#include "AsciiDataWriter.hh"
#include "Exceptions.hh"

#include <string.h>
#include <mpi.h>

AsciiDataWriter::AsciiDataWriter(int rank, int size)
  : DataWriter(rank, size), time_dim_(true)
{}

AsciiDataWriter::AsciiDataWriter(int rank, int size, 
                                 const char *fn, Result &result)
  : DataWriter(rank, size), time_dim_(true)
{
  set_filename(fn);
  add_variable(result);
}

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
    throw DataWriterException("Result must have at least one dimension.");
  
  dim_len_ = 0;
  for (unsigned int i = 0; i < dim_lens.size(); i++)
    dim_len_ += dim_lens[i];

  time_dim_ = result.has_time_dimension();
}

unsigned int AsciiDataWriter::write_data(unsigned int time_step, 
                                         Data &data, MPI_Datatype t, 
                                         void *ptr, unsigned int len)
{
  if (!file_.is_open()) 
  {
    init();
    
    if (!file_.is_open()) 
      throw DataWriterException("Unable to open output file!");
  }

  // If this result doesn't have a time dimension, then move to the
  // start of the file before writing anything; overwrite existing
  // data.
  if (!time_dim_)
    file_.seekp(0);

  unsigned int tbw = 0, bytes_written = 0;

  while (tbw < len) {
    bytes_written = write_data(t, ptr, len);
    tbw += bytes_written;

    // Kind of ugly, but necessary to work on AIX
    char *byte_ptr = static_cast<char *>(ptr);
    byte_ptr += bytes_written;
    ptr = static_cast<void *>(byte_ptr);

    if (!time_dim_)
      file_ << endl;
  }

  file_ << endl;

  // TEMPORARY! PERFORMANCE KILL (probably)
  file_.flush();

  return tbw;
}

unsigned int AsciiDataWriter::write_data(MPI_Datatype t, void *ptr, 
                                  unsigned int len)
{
  unsigned int bytes_written = 0; 

  if (t == MPI_CHAR) 
    bytes_written = write_array<signed char>(ptr, len);
  else if (t == MPI_SHORT)
    bytes_written = write_array<signed short int>(ptr, len);
  else if (t == MPI_INT)
    bytes_written = write_array<signed int>(ptr, len);
  else if (t == MPI_LONG)
    bytes_written = write_array<signed long int>(ptr, len);
  else if (t == MPI_UNSIGNED_CHAR)
    bytes_written = write_array<unsigned char>(ptr, len);
  else if (t == MPI_UNSIGNED_SHORT)
    bytes_written = write_array<unsigned short int>(ptr, len);
  else if (t == MPI_UNSIGNED)
    bytes_written = write_array<unsigned int>(ptr, len);
  else if (t == MPI_UNSIGNED_LONG)
    bytes_written = write_array<unsigned long int>(ptr, len);
  else if (t == MPI_FLOAT)
    bytes_written = write_array<float>(ptr, len);
  else if (t == MPI_DOUBLE)
    bytes_written = write_array<double>(ptr, len);
  else if (t == MPI_LONG_DOUBLE)
    bytes_written = write_array<long double>(ptr, len);
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
      // Get the size of the data type, so we can send the right
      // number for the third arg of write_data. 
      int sz;
      MPI_Type_size(dts[i], &sz);
      bytes_written = write_data(dts[i], ptr, ints[i] * sz);

      // Kind of ugly, but necessary to work on AIX
      char *byte_ptr = static_cast<char *>(ptr);
      byte_ptr += bytes_written;
      ptr = static_cast<void *>(byte_ptr);
      i++;
    }

    delete[] ints;
    delete[] aints;
    delete[] dts;
  }
  
  return bytes_written;
}

template <class T>
unsigned int AsciiDataWriter::write_array(void *ptr, unsigned int len)
{
  unsigned int items_avail = items_avail = len / sizeof(T);
  T *r = static_cast<T *>(ptr);

  for (unsigned int i = 0; i < items_avail; i++)
  {
    file_ << *(r++) << "\t";
  }

  return len;
}
