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

#include "MatlabDataWriter.hh"
#include "Exceptions.hh"

#include <string.h>
#include <mpi.h>
#include <stdio.h>
#include <time.h>

unsigned int sizeof_matlab(MATLAB_data_type type)
{
  switch (type)
  {
  case miINT8:
  case miUINT8:
    return 1;
    break;

  case miINT16:
  case miUINT16:
    return 2;
    break;

  case miINT32:
  case miUINT32:
    return 4;
    break;

  case miINT64:
  case miUINT64:
    return 8;
    break;

  case miSINGLE:
    return sizeof(float);
    break;

  case miDOUBLE:
    return sizeof(double);
    break;
  }

  return 0;
}

MATLAB_data_type MPI_to_matlab_dt(MPI_Datatype mpi_type)
{
  if (mpi_type == MPI_CHAR)
    return miINT8;

  if (mpi_type == MPI_BYTE)
    return miUINT8;

  if (mpi_type == MPI_SHORT)
    return miINT16;

  if (mpi_type == MPI_INT)
    return miINT32;

  if (mpi_type == MPI_LONG)
    return miINT64;

  if (mpi_type == MPI_FLOAT)
    return miSINGLE;

  if (mpi_type == MPI_DOUBLE)
    return miDOUBLE;

  if (mpi_type == MPI_UNSIGNED_CHAR)
    return miUINT8;

  if (mpi_type == MPI_UNSIGNED_SHORT)
    return miUINT16;

  if (mpi_type == MPI_UNSIGNED)
    return miUINT32;

  if (mpi_type == MPI_UNSIGNED_LONG)
    return miUINT64;

  return miINT8;
}

/************************************************************
 * MatlabElement function implementations
 ************************************************************/
MatlabElement::MatlabElement()
  : buffer_(0), buffer_pos_(0), buffer_size_(0)
{
  tag_.datatype = miINT8;
  tag_.num_bytes = 0;
}

MatlabElement::MatlabElement(MATLAB_data_type type)
  : buffer_(0), buffer_pos_(0), buffer_size_(0)
{
  tag_.datatype = type;
  tag_.num_bytes = 0;
}

MatlabElement::~MatlabElement()
{
  if (buffer_)
    delete[] buffer_;
}

unsigned int MatlabElement::get_num_bytes()
{
  unsigned int padding = (8 - (tag_.num_bytes) % 8);
  if (padding >= 8)
    padding = 0;

  return tag_.num_bytes + 8 + padding;
}

// void MatlabElement::write_data(ostream &stream, unsigned int num_bytes,
//                                const void *ptr)
// {

// }

void MatlabElement::append_buffer(unsigned int num_bytes, const void *ptr)
{
  if (!buffer_)
  {
    buffer_ = new char[1024];
    buffer_size_ = 1024;
    
    if (!buffer_)
      throw MemoryException();
  }

  if (buffer_pos_ + num_bytes > buffer_size_)
  {
    buffer_size_ = buffer_pos_ + num_bytes + 1024;
    char *new_buf = new char[buffer_size_];
    if (!new_buf)
      throw MemoryException();

    memmove(new_buf, buffer_, tag_.num_bytes);
    delete[] buffer_;
    buffer_ = new_buf;
  }

  memmove(buffer_ + buffer_pos_, static_cast<const char *>(ptr), num_bytes);
  buffer_pos_ += num_bytes;
  tag_.num_bytes += num_bytes;
}

void MatlabElement::write_buffer(ostream &stream)
{
  stream.write(reinterpret_cast<char *>(&tag_), 8);

  // Rejigger the data so that it is in column major order, otherwise
  // the matrix in matlab will have to be transposed.
  stream.write(buffer_, tag_.num_bytes);

  int padding = 8 - (tag_.num_bytes) % 8;
  if (padding > 0 && padding < 8)
    for (int i = 0; i < padding; i++)
      stream.put(0);
  stream.flush();
}

// void MatlabElement::update_file_offset(int bytes)
// {

// }

/************************************************************
 * MatlabArray function implementations
 ************************************************************/
MatlabArray::MatlabArray(const char *name, 
                         const vector<int> &dim_lens, bool time_dim, 
                         MATLAB_data_type type,
                         bool complex)
{
  name_ = name;
  memset(static_cast<void *>(&flags_), 0, sizeof(array_flags_t));
  flags_.flags.bits.complex_bit = complex ? 1 : 0;
  flags_.flags.bits.global_bit = 0; // 1;
  flags_.flags.bits.logical_bit = 0;
  //flags_.flags.byte = 0x8 + 0x4;
  flags_.array_class = static_cast<uint8_t>(get_array_class(type));

  time_dim_ = time_dim;
    
  // Array flags
  me_flags_.set_type(miUINT32);
  me_flags_.append_buffer(8, reinterpret_cast<void *>(&flags_));

  // Dimension lengths
  me_dim_lens_.set_type(miINT32);
  num_dims_ = dim_lens.size();

  if (time_dim)
    num_dims_++;

  dim_lengths_ = new int32_t[num_dims_];
  
  vector<int>::const_iterator iter;
  vector<int>::const_iterator iter_e = dim_lens.end();
  int i = 0;
  int sz = 0;
  
  if (time_dim)
  {
    i++;
    dim_lengths_[0] = 0;
  }

  for(iter = dim_lens.begin(); iter != iter_e; ++iter, ++i)
  {
    sz += *iter;
    dim_lengths_[i] = *iter;
  }
  
  // Array name
  me_name_.set_type(miINT8);
  me_name_.append_buffer(name_.length(), 
                         static_cast<const void *>(name));

  // Real data
  me_data_.set_type(type);

  // Optional Imaginary data ... 

  tag_.datatype = miMATRIX;
}

MatlabArray::~MatlabArray()
{
  if (dim_lengths_)
    delete[] dim_lengths_;
}

void MatlabArray::append_buffer(unsigned int num_bytes, const void *ptr)
{
  me_data_.append_buffer(num_bytes, ptr);
  tag_.num_bytes += num_bytes;

  if (time_dim_)
    dim_lengths_[0] += 1;
}

void MatlabArray::write_buffer(ostream &stream)
{
  me_dim_lens_.append_buffer(sizeof(int32_t) * num_dims_,
                             static_cast<const void *>(dim_lengths_));

  tag_.num_bytes = 0;
  tag_.num_bytes += me_data_.get_num_bytes();
  tag_.num_bytes += me_name_.tag_.num_bytes + 4;
  int padding = 8 - (me_name_.tag_.num_bytes + 4) % 8;

  if (padding == 8)
    padding = 0;

  tag_.num_bytes += padding;
  tag_.num_bytes += me_flags_.get_num_bytes();
  tag_.num_bytes += me_dim_lens_.get_num_bytes();

  stream.write(reinterpret_cast<char *>(&tag_), sizeof(data_tag_t));
  stream.flush();

  me_flags_.write_buffer(stream);
  me_dim_lens_.write_buffer(stream);

  // The writing of the name doesn't work the same way as for other elements. 
  uint16_t temp = static_cast<uint16_t>(me_name_.tag_.num_bytes);
  stream.write(reinterpret_cast<char *>(&temp), 2);
  stream.flush();
  
  temp = static_cast<uint16_t>(me_name_.tag_.datatype);
  stream.write(reinterpret_cast<char *>(&temp), 2);
  stream.flush();
  
  stream.write(me_name_.buffer_, me_name_.tag_.num_bytes);
  stream.flush();

  if (padding > 0 && padding < 8)
    for (int i = 0; i < padding; i++)
      stream.put(0);
  stream.flush();

  me_data_.write_buffer(stream);
}

// void MatlabArray::write_data(ostream &stream, unsigned int num_bytes,
//                              const void *ptr, bool buffer)
// {
  
  
// }

unsigned int MatlabArray::get_num_bytes()
{
  return tag_.num_bytes;
}

MATLAB_array_type MatlabArray::get_array_class(MATLAB_data_type type)
{
  MATLAB_array_type ret = mxDOUBLE_CLASS;

  if (type >= miINT8 && type < miSINGLE)
    ret = static_cast<MATLAB_array_type>(type + mxSINGLE_CLASS);
  else if (type == miSINGLE)
    ret = static_cast<MATLAB_array_type>(mxSINGLE_CLASS);
  else if (type == miDOUBLE)
    ret = static_cast<MATLAB_array_type>(mxDOUBLE_CLASS);
  else
    throw DataWriterException("Invalid MATLAB data type given.");

  return ret;
}

/************************************************************
 * MatlabDataWriter function implementations
 ************************************************************/

MatlabDataWriter::MatlabDataWriter(int rank, int size)
  : DataWriter(rank, size), buffered_(true)
{
  header_setup();
}

MatlabDataWriter::MatlabDataWriter(int rank, int size, 
                                 const char *fn, Result &result)
  : DataWriter(rank, size), buffered_(true)
{
  header_setup();
  filename_ = fn;
  add_variable(result);
}

MatlabDataWriter::~MatlabDataWriter()
{
}

void MatlabDataWriter::header_setup()
{
  time_t t = time(0);
  memset(static_cast<void *>(&header_), ' ', 124);
  snprintf(header_.text, 124, "MATLAB 5.0 MAT-file, Platform: %s, Created on: %s by phred %s", 
           PLATFORM, ctime(&t), PACKAGE_VERSION);
  header_.version = 0x0100;
  //header_.version = 0x0001;
  header_.endian.bytes[0] = 'M';
  header_.endian.bytes[1] = 'I';
}

void MatlabDataWriter::init(const Grid &grid)
{
  //test();

  if (filename_.length() > 0 && rank_ == 0)
  {
    file_.open(filename_.c_str(), ofstream::out | ofstream::binary
               | ofstream::trunc);
  }
}

void MatlabDataWriter::deinit(const Grid &grid)
{
  // If data is buffered, and things are ready to write, write them. 
  map<string, MatlabArray *>::iterator iter;
  map<string, MatlabArray *>::iterator iter_e = vars_.end();

  if (rank_ == 0 && file_.is_open())
  {
    file_.write((char *)(&header_), 128);

    for(iter = vars_.begin(); iter != iter_e; ++iter)
      iter->second->write_buffer(file_);

    file_.close();
  }

  for (iter = vars_.begin(); iter != iter_e; ++iter)
  {
    delete iter->second;
  }
  
  vars_.clear();
}

void MatlabDataWriter::add_variable(Result &result)
{
  const map<string, Variable *> vars = result.get_variables();
  map<string, Variable *>::const_iterator iter;
  map<string, Variable *>::const_iterator iter_e = vars.end();
  
  for (iter = vars.begin(); iter != iter_e; ++iter)
  {
    Variable *var = iter->second;
    const vector<int> &dim_lens = var->get_dim_lengths();
    
    if (dim_lens.size() == 0)
      throw DataWriterException("Result must have at least one dimension.");
    
    vars_[var->get_name()] = 
      new MatlabArray(var->get_name().c_str(), 
                      dim_lens, var->has_time_dimension(), 
                      MPI_to_matlab_dt(var->get_element_type()),
                                       false);
  }
}

unsigned int MatlabDataWriter::write_data(unsigned int time_step, 
                                         Variable &variable, MPI_Datatype t, 
                                         void *ptr, unsigned int len)
{
  const Data &data = variable.get_data();
    
  if (!file_.is_open()) 
    throw DataWriterException("MatlabDataWriter: File should already be open!");

  try {
    vars_[variable.get_name()]->append_buffer(len, ptr);
  } catch (MemoryException e) {
    if (rank_ == 0 && file_.is_open())
    {
      file_.write(reinterpret_cast<char *>(&header_), 128);
      
      map<string, MatlabArray *>::iterator iter;
      map<string, MatlabArray *>::iterator iter_e = vars_.end();

      for(iter = vars_.begin(); iter != iter_e; ++iter)
        iter->second->write_buffer(file_);
      
      file_.close();
    }
    throw MemoryException();
  }
  return len;
}

void MatlabDataWriter::test()
{
  vector<int> dims, dims2;
  dims.push_back(2);
  dims.push_back(2);
  dims2.push_back(4);
  dims2.push_back(2);
  dims2.push_back(2);
  double data2[] = {1, 2, 3, 4, 1, 2, 3, 4, 7, 8, 9, 5, 7, 8, 9, 5};

  MatlabArray *ma = new MatlabArray("ab", dims, false, miDOUBLE, false);
  MatlabArray *ma2 = new MatlabArray("test", dims2, false, miDOUBLE, false);

  double data[] = {1, 2, 3, 4};
  ma->append_buffer(4 * sizeof(double), reinterpret_cast<void *>(data));
  ma2->append_buffer(16 * sizeof(double), reinterpret_cast<void *>(data2));
  
  ofstream tf("a.mat",  ofstream::out | ofstream::binary
             | ofstream::trunc);

  tf.write(reinterpret_cast<char *>(&header_), 128);
  tf.flush();
  ma2->write_buffer(tf);
  delete ma2;

  ma->write_buffer(tf);
  delete ma;

  tf.close();
}
