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

#include "config.h"

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
MatlabElement::MatlabElement(bool compress)
  : buffer_(0), buffer_pos_(0), buffer_size_(0), compress_(compress)
{
  tag_.datatype = miINT8;
  tag_.num_bytes = 0;
}

MatlabElement::MatlabElement(MATLAB_data_type type, bool compress)
  : buffer_(0), buffer_pos_(0), buffer_size_(0), compress_(compress)
{
  tag_.datatype = type;
  tag_.num_bytes = 0;
}  

MatlabElement::~MatlabElement()
{
  if (buffer_)
    delete[] buffer_;
}

void MatlabElement::set_type(MATLAB_data_type type)
{
  source_datatype_ = type;
  tag_.datatype = type;
  tag_.num_bytes = 0;
}

void MatlabElement::compress()
{
  if (compress_)
  {
    if (tag_.datatype <= 9)
    {
      switch (source_datatype_)
      {
      case miINT8:
        compress_helper<signed char>();
        break;
      case miUINT8:
        compress_helper<unsigned char>();
        break;
      case miINT16:
        compress_helper<signed short>();
        break;
      case miUINT16:
        compress_helper<unsigned short>();
        break;
      case miINT32:
        compress_helper<signed int>();
        break;
      case miUINT32:
        compress_helper<unsigned int>();
        break;
      case miSINGLE:
        compress_helper<float>();
        break;
      case miDOUBLE:
        compress_helper<double>();
        break;
      }
    } else
      compressed_tag_ = tag_;
  }
}

void MatlabElement::write_compress(ostream &stream)
{
  if (compress_)
  {
    if (tag_.datatype <= 9)
    {
      switch (source_datatype_)
      {
      case miINT8:
        write_compress_helper<signed char>(stream);
        break;
      case miUINT8:
        write_compress_helper<unsigned char>(stream);
        break;
      case miINT16:
        write_compress_helper<signed short>(stream);
        break;
      case miUINT16:
        write_compress_helper<unsigned short>(stream);
        break;
      case miINT32:
        write_compress_helper<signed int>(stream);
        break;
      case miUINT32:
        write_compress_helper<unsigned int>(stream);
        break;
      case miSINGLE:
        write_compress_helper<float>(stream);
        break;
      case miDOUBLE:
        write_compress_helper<double>(stream);
        break;
      }
    } else {
    throw DataWriterException("MatlabElement::write_compress: can't write this data type!");
    }

  } else
    throw DataWriterException("MatlabElement::write_compress erroniously called!");
}

unsigned int MatlabElement::get_num_bytes()
{
  unsigned int nbytes = 0;
  data_tag_t &tag = tag_;

  if (compress_)
  {
    compress();
    tag = compressed_tag_;
  }

  if (tag.num_bytes <= 4)
  {
    nbytes = 8;
  } else {
    unsigned int padding = (8 - (tag.num_bytes) % 8);
    if (padding >= 8)
      padding = 0;

    nbytes = tag.num_bytes + 8 + padding;
  }

  return nbytes;
}

void MatlabElement::overwrite_buffer(unsigned int num_bytes, const void *ptr)
{
  if (!buffer_) 
  {
    buffer_size_ = num_bytes;
    buffer_ = new char[buffer_size_];
    
    if (!buffer_)
      throw MemoryException();

    memset(buffer_, 0, sizeof(char) * buffer_size_);
  }

  if (num_bytes != buffer_size_)
    throw DataWriterException("MatlabDataWriter: non-time variable buffer size should not change!");

  memmove(buffer_, ptr, buffer_size_);

  buffer_pos_ = num_bytes;
  tag_.num_bytes = num_bytes;
}

void MatlabElement::append_buffer(unsigned int num_bytes, const void *ptr)
{
  if (!buffer_) 
  {
    buffer_size_ = num_bytes;
    buffer_ = new char[buffer_size_];
    
    if (!buffer_)
      throw MemoryException();

    memset(buffer_, 0, sizeof(char) * buffer_size_);
  }

  if (buffer_pos_ + num_bytes > buffer_size_)
  {
    buffer_size_ = buffer_pos_ + num_bytes + 1024;
    char *new_buf = new char[buffer_size_];
    if (!new_buf)
      throw MemoryException();
    memset(new_buf, 0, sizeof(char) * buffer_size_);

    memmove(new_buf, buffer_, tag_.num_bytes);
    delete[] buffer_;
    buffer_ = new_buf;
  }

  memmove(buffer_ + buffer_pos_, ptr, num_bytes);
  buffer_pos_ += num_bytes;
  tag_.num_bytes += num_bytes;
}

void MatlabElement::write_buffer(ostream &stream)
{
  data_tag_t tag = tag_;
  int padding = 0;

  if (compress_)
    tag = compressed_tag_;

  if (tag.num_bytes <= 4)
  {
#ifdef WORDS_BIGENDIAN
    uint16_t temp = static_cast<uint16_t>(tag.num_bytes);
    stream.write(reinterpret_cast<char *>(&temp), 2);
    
    temp = static_cast<uint16_t>(tag.datatype);
    stream.write(reinterpret_cast<char *>(&temp), 2);
#else
    uint16_t temp = static_cast<uint16_t>(tag.datatype);
    stream.write(reinterpret_cast<char *>(&temp), 2);
    
    temp = static_cast<uint16_t>(tag.num_bytes);
    stream.write(reinterpret_cast<char *>(&temp), 2);    
#endif
 
    padding = 8 - (tag.num_bytes + 4) % 8;

  } else {

    stream.write(reinterpret_cast<char *>(&tag), 8);
    padding = 8 - (tag.num_bytes) % 8;
  }

  if (compress_)
    write_compress(stream);
  else
    stream.write(buffer_, tag.num_bytes);
    
  if (padding > 0 && padding < 8)
    for (int i = 0; i < padding; i++)
      stream.put(0);
  
  stream.flush();
}

void MatlabElement::reshape_buffer(int N, int M, MPI_Datatype type)
{
  MPI_Datatype new_t, flat_t, new2_t;
  MPI_Status status;

  char *new_buf = new char[buffer_size_];
  
#ifdef DEBUG
  //cerr << "Reshaping buffer to " << N << " by " << M << endl;
#endif

  if (!new_buf)
    throw MemoryException();

  memset(new_buf, 0, sizeof(char) * buffer_size_);

  int type_size;
  MPI_Type_size(type, &type_size);

  MPI_Type_vector(N, 1, M, type, &new_t);
  MPI_Type_commit(&new_t);

  MPI_Type_hvector(M, 1, type_size, new_t, &new2_t);
  MPI_Type_commit(&new2_t);

  MPI_Type_contiguous(N * M, type, &flat_t);
  MPI_Type_commit(&flat_t);

  MPI_Sendrecv(buffer_, 1, new2_t, 0, 0, 
               new_buf, 1, flat_t, 0, 0, 
               MPI_COMM_WORLD, &status);

  delete[] buffer_;
  buffer_ = new_buf;

  MPI_Type_free(&new_t);
  MPI_Type_free(&new2_t);
  MPI_Type_free(&flat_t);
}

/************************************************************
 * MatlabArray function implementations
 ************************************************************/
MatlabArray::MatlabArray(const char *name, 
                         const vector<int> &dim_lens, bool time_dim, 
                         MPI_Datatype type,
                         bool complex)
{
  vector<int>::const_iterator iter;
  vector<int>::const_iterator iter_e = dim_lens.end();
  int i = 0;
  int sz = 0;
  
  name_ = name;
  time_dim_ = time_dim;
  mpi_type_ = type;
    
  memset(reinterpret_cast<void *>(&flags_), 0, sizeof(array_flags_t));

#ifdef WORDS_BIGENDIAN
  flags_.flags.flag_bytes[3] = 
    static_cast<uint8_t>(get_array_class(MPI_to_matlab_dt(type)));

  flags_.flags.flag_bytes[2] = complex ? 0x8 : 0x0;
  //flags_.flags.flag_bytes[2] = flags_.flags.flag_bytes[2] | 0x4; // global
  //flags_.flags.flag_bytes[2] = flags_.flags.flag_bytes[2] | 0x2; // logical 
#else
  flags_.flags.flag_bytes[0] = 
    static_cast<uint8_t>(get_array_class(MPI_to_matlab_dt(type)));

  flags_.flags.flag_bytes[1] = complex ? 0x8 : 0x0;
  //flags_.flags.flag_bytes[1] = flags_.flags.flag_bytes[2] | 0x4; // global
  //flags_.flags.flag_bytes[1] = flags_.flags.flag_bytes[2] | 0x2; // logical 
#endif

  // Array flags
  me_flags_.set_type(miUINT32);
  me_flags_.append_buffer(8, reinterpret_cast<void *>(&flags_));

  // Dimension lengths
  me_dim_lens_.set_type(miINT32);
  num_dims_ = dim_lens.size();

  if (time_dim)
    num_dims_++;

  if (num_dims_ < 2)
  {
    i++;
    num_dims_++;
  }

  dim_lengths_ = new int32_t[num_dims_];

  if (!dim_lengths_)
    throw MemoryException();

  if (i > 0)
    dim_lengths_[0] = 1;
  
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

  me_dim_lens_.append_buffer(sizeof(int32_t) * num_dims_,
                             static_cast<const void *>(dim_lengths_));

  // Array name
  me_name_.set_type(miINT8);
  me_name_.append_buffer(name_.length(), 
                         static_cast<const void *>(name));

  // Real data
  me_data_.set_type(MPI_to_matlab_dt(type));

  // Optional Imaginary data ... 


  // Setup our tag...
  tag_.datatype = miMATRIX;
  tag_.num_bytes = 0;
  tag_.num_bytes += me_name_.get_num_bytes();
  tag_.num_bytes += me_data_.get_num_bytes();
  tag_.num_bytes += me_flags_.get_num_bytes();
  tag_.num_bytes += me_dim_lens_.get_num_bytes();

}

MatlabArray::~MatlabArray()
{
  if (dim_lengths_)
    delete[] dim_lengths_;
}

void MatlabArray::append_buffer(unsigned int num_bytes, const void *ptr)
{
  if (time_dim_) 
  {
    me_data_.append_buffer(num_bytes, ptr);
    tag_.num_bytes += num_bytes;

    dim_lengths_[0] += 1;

    me_dim_lens_.overwrite_buffer(sizeof(int32_t) * num_dims_,
                                  static_cast<const void *>(dim_lengths_));

  } else {
    me_data_.overwrite_buffer(num_bytes, ptr);
  }
}

void MatlabArray::write_buffer(ostream &stream)
{
  // Find the smallest data type that can be used to store the data on
  // disk. This is matlab's "compressed" format. 

  get_num_bytes();
  stream.write(reinterpret_cast<char *>(&tag_), sizeof(data_tag_t));
  stream.flush();

  me_flags_.write_buffer(stream);
  me_dim_lens_.write_buffer(stream);
  me_name_.write_buffer(stream);

  if (time_dim_)
  {
    unsigned int cols = 0;
    for (int i = 1; i < num_dims_; i++)
      cols += dim_lengths_[i];

    me_data_.reshape_buffer(dim_lengths_[0], cols, mpi_type_);
  }

  me_data_.write_buffer(stream);
}

unsigned int MatlabArray::get_num_bytes()
{
  tag_.num_bytes = 0;
  tag_.num_bytes += me_name_.get_num_bytes();
  tag_.num_bytes += me_data_.get_num_bytes();
  tag_.num_bytes += me_flags_.get_num_bytes();
  tag_.num_bytes += me_dim_lens_.get_num_bytes();

  return tag_.num_bytes;
}

MATLAB_array_type MatlabArray::get_array_class(MATLAB_data_type type)
{
  MATLAB_array_type ret = mxDOUBLE_CLASS;

  if (type >= miINT8 && type < miSINGLE)
    ret = static_cast<MATLAB_array_type>(type + mxSINGLE_CLASS);
  else if (type == miSINGLE)
    ret = static_cast<MATLAB_array_type>(mxSINGLE_CLASS);
  //ret = static_cast<MATLAB_array_type>(mxDOUBLE_CLASS);
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

#ifdef WORDS_BIGENDIAN
  header_.endian[0] = 'M';
  header_.endian[1] = 'I';
#else
  header_.endian[1] = 'M';
  header_.endian[0] = 'I';
#endif
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
    const vector<Dimension> &dimensions = var->get_dimensions();
    
    if (dimensions.size() == 0)
      throw DataWriterException("Result must have at least one dimension.");
    
    vector<int> dim_lens;
    for(unsigned int idx; idx < dimensions.size(); idx++)
      dim_lens.push_back(dimensions[idx].global_len_);

    vars_[var->get_name()] = 
      new MatlabArray(var->get_name().c_str(), 
                      dim_lens, var->has_time_dimension(), 
                      var->get_element_type(),
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
    if (len > 0)
    {
#ifdef DEBUG
//       cerr << "Buffering data for matlab output from " << 
//         variable.get_name() << ", " << len << " bytes. " << endl;
//       cerr << "\tMatlabDataWriter, recieved data, floats are: " 
//            << reinterpret_cast<float *>(ptr)[0] << " and " 
//            << reinterpret_cast<float *>(ptr)[1] << endl;
//       cerr << "\tPointer is " << ptr << endl;
#endif

      vars_[variable.get_name()]->append_buffer(len, ptr);
    }
  } catch (MemoryException e) {
    if (rank_ == 0 && file_.is_open())
    {
      file_.write(reinterpret_cast<char *>(&header_), 128);
      
      map<string, MatlabArray *>::iterator iter;
      map<string, MatlabArray *>::iterator iter_e = vars_.end();

      for(iter = vars_.begin(); iter != iter_e; ++iter)
      {
        iter->second->write_buffer(file_);
      }

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
  dims2.push_back(4);
  dims2.push_back(2);
  dims2.push_back(2);
  float data2[] = {1, 2, 3, 4, 1, 2, 3, 4, 7, 8, 9, 5, 7, 8, 9, 5};
  short data3[] = {1, 2, 3, 12, 43, 45, 87, 102, 123};

  MatlabArray *ma = new MatlabArray("test2", dims, true, MPI_DOUBLE, false);
  MatlabArray *ma2 = new MatlabArray("abc", dims2, false, MPI_FLOAT, false);
  MatlabArray *ma3 = new MatlabArray("shortdata", dims, true, MPI_SHORT, false);
  if (!ma || !ma3 || !ma2)
    throw MemoryException();

  double data[] = {2, 3, 4, 5};
  ma->append_buffer(2 * sizeof(double), reinterpret_cast<void *>(data));
  ma->append_buffer(2 * sizeof(double), reinterpret_cast<void *>(data + 2));
  ma->append_buffer(2 * sizeof(double), reinterpret_cast<void *>(data));

  ma2->append_buffer(16 * sizeof(float), reinterpret_cast<void *>(data2));
  ma2->append_buffer(16 * sizeof(float), reinterpret_cast<void *>(data2));
  
  ma3->append_buffer(2 *sizeof(short), reinterpret_cast<void *>(data3));
  ma3->append_buffer(2 *sizeof(short), reinterpret_cast<void *>(data3 + 2));
  ma3->append_buffer(2 *sizeof(short), reinterpret_cast<void *>(data3));

  ofstream tf("a.mat",  ofstream::out | ofstream::binary
             | ofstream::trunc);

  tf.write(reinterpret_cast<char *>(&header_), 128);
  tf.flush();

  ma2->write_buffer(tf);
  delete ma2;

  ma->write_buffer(tf);
  delete ma;

  ma3->write_buffer(tf);
  delete ma3;

  tf.close();
}
