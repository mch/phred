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

#ifndef MATLAB_DATA_WRITER_H
#define MATLAB_DATA_WRITER_H

#include "DataWriter.hh"
#include <fstream>
#include <exception>
#include <map>
#include <string>

#include <sys/types.h>

#ifdef __GNUC__
#include <stdint.h> // uint32_t on linux
#endif

#include <mpi.h>
#include <limits.h>

using namespace std;

/**
 * MATLAB output data types
 */
enum MATLAB_data_type {
  miINT8 = 1, /**< 8 bit signed */
  miUINT8 = 2, /**< 8 bit unsigned */
  miINT16 = 3, /**< 16 bit signed */
  miUINT16 = 4, /**< 16 bit unsigned */
  miINT32 = 5, /**< 32 bit signed */
  miUINT32 = 6, /**< 32 bit unsigned */
  miSINGLE = 7, /**< IEEE 754 single format */
  miDOUBLE = 9, /**< IEEE 754 double format */
  miINT64 = 12, /**< 64 bit signed */
  miUINT64 = 13, /**< 64 bit unsigned */
  miMATRIX = 14, /**< Matrix */ 
};

/**
 * A little helper to convert MPI data types to MATLAB data types. 
 *
 * @param mpi_type The MPI Data type to convert
 * @return MATLAB_data_type
 */ 
MATLAB_data_type MPI_to_matlab_dt(MPI_Datatype mpi_type);

/**
 * Tells the number of bytes for a given MATLAB data type.
 */ 
unsigned int sizeof_matlab(MATLAB_data_type type);

/**
 * MATLAB Array classes 
 */ 
enum MATLAB_array_type {
  mxCELL_CLASS = 1, /**< Cell array */
  mxSTRUCT_CLASS = 2, /**< Structure */
  mxOBJECT_CLASS = 3, /**< Object */
  mxCHAR_CLASS = 4, /**< Character array */
  mxSPARSE_CLASS = 5, /**< Sparse array */
  mxDOUBLE_CLASS = 6,  /**< Double precision array */
  mxSINGLE_CLASS = 7, /**< Single precision array */
  mxINT8_CLASS = 8, /**< 8 bit signed integer */
  mxUINT8_CLASS = 9, /**< 8 bit unsigned integer */
  mxINT16_CLASS = 10, /**< 16 bit signed integer */
  mxUINT16_CLASS = 11, /**< 16 bit unsigned integer */
  mxINT32_CLASS = 12, /**< 32 bit signed integer */
  mxUINT32_CLASS = 13, /**< 32 bit unsigned integer */
};

/**
 * Data tag
 */
typedef struct {
  uint32_t datatype;
  uint32_t num_bytes;
} data_tag_t;    

/**
 * Array flags format
 */
typedef struct {
  union flags_t {
    uint32_t flags32;
    uint8_t flag_bytes[4];
  } flags;
  uint32_t undefined;
} array_flags_t;

/**
 * Represents a data element in a MATLAB file. 
 *
 * \bug This buffers data until we can write it out. The buffer
 * implementation is extremly poor; it reallocates a new chunk of
 * memory when it is needed and just copies everything. Ugg. Chaining
 * would probably be better.
 */
class MatlabElement {
  friend class MatlabArray;
  friend class MatlabDataWriter;
private:
  MatlabElement(MATLAB_data_type type, bool compress = false);
  MatlabElement(bool compress = false);

  // Not allowed. 
  MatlabElement(const MatlabElement &rhs)
  {}
  const MatlabElement &operator=(const MatlabElement &rhs)
  {
    return *this;
  }

protected:
  data_tag_t tag_; /**< Type and size of data */ 
  data_tag_t compressed_tag_; /**< Type and size of data if written in
                                 compressed format */

  unsigned char num_pad_bytes_; /**< Number of pad bytes on the end of
                                   the data to make it align to 64 bit
                                   boundaries. */  

  MATLAB_data_type source_datatype_; /**< Indicates the matlab
                                        datatype of the buffered
                                        data. */

  unsigned int file_offset_; /**< Position (in bytes) in the file of
                                the start of the tag for this element */ 

  char *buffer_; /**< Buffer for data, so we can just write it all out
                    at the end. */ 
  unsigned int buffer_pos_; /**< Place to start appending data in the
                               buffer. */ 
  unsigned int buffer_size_; /**< Actual allocated memory */ 

  bool compress_; /**< Report byte counts and write data in MATLAB's
                     'compressed' mode, which simply uses the smallest possible
                     data-type to represent output. */
  /**
   * Check the data and compute the size and type of the data as it
   * would be written. If compress_ is false, this function does
   * nothing. Otherwise, it finds the smallest data type that may be
   * used to represent the buffered data, and stores that in compressed_tag_. 
   */
  void compress();

  template<class T>
  void compress_helper()
  {
    T *ptr = static_cast<T *>(buffer_);
    unsigned int size = tag_.num_bytes / sizeof(T);
    bool is_int = true;
    double max = 0;
    double min = 0;
    
    for (unsigned int idx = 0; idx < size; idx++)
    {      
      T temp = ptr[idx];
      if (temp > max)
        max = temp;
      if (temp < min)
        min = temp;
      if (temp % 1 > 0)
        is_int = false;
    }

    if (is_int) {
      if (min >= UCHAR_MIN && max <= UCHAR_MAX)
      {
        compressed_tag_.datatype = miUINT8;
        compressed_tag_.num_bytes = size * sizeof(unsigned char);
      }
      else if (min >= CHAR_MIN && max <= CHAR_MAX)
      {
        compressed_tag_.datatype = miINT8;
        compressed_tag_.num_bytes = size * sizeof(signed char);
      }
      else if (min >= 0 && max <= USHRT_MAX)
      {
        compressed_tag_.datatype = miUINT16;
        compressed_tag_.num_bytes = size * sizeof(unsigned short);
      }
      else if (min >= SHRT_MIN && max <= SHRT_MAX)
      {
        compressed_tag_.datatype = miINT16;
        compressed_tag_.num_bytes = size * sizeof(signed short);
      }
      else if (min >= 0 && max <= UINT_MAX)
      {
        compressed_tag_.datatype = miUINT32;
        compressed_tag_.num_bytes = size * sizeof(unsigned int);
      }
      else if (min >= INT_MIN && max <= INT_MAX)
      {
        compressed_tag_.datatype = miINT32;
        compressed_tag_.num_bytes = size * sizeof(signed int);
      }
    } 
    else
    {
      compressed_tag_ = tag_;
    }
  }

public:
  virtual ~MatlabElement();

  /**
   * How many bytes (data and tag) will this variable write? 
   */
  virtual unsigned int get_num_bytes();

  /**
   * Set the MATLAB data type for this element. 
   *
   * @param type The MATLAB datatype
   */
  void set_type(MATLAB_data_type type);

  /**
   * Overwrite data in the buffer
   *
   * @param num_bytes The number of bytes to write
   * @param ptr The memory to get the data from
   */ 
  virtual void overwrite_buffer(unsigned int num_bytes, 
                                const void *ptr);

  /**
   * Append to the buffer. 
   *
   * @param num_bytes The number of bytes to append to the end of
   * this element. 
   * @param ptr The memory to get the data from. 
   */
  virtual void append_buffer(unsigned int num_bytes, 
                             const void *ptr);

  /**
   * Write the buffered data to disk, along with it's tag telling the
   * data type and length.
   *
   * @param stream the stream to write the data to
   */ 
  virtual void write_buffer(ostream &stream);

  /**
   * Reshape the data into N rows and M columns.
   * 
   * @param N number of rows
   * @param M Number of columns
   */
  virtual void reshape_buffer(int N, int M, MPI_Datatype type);
  
};

/**
 * Represents a chunk of array data, like a matrix, in a MATLAB
 * file. Knows how to write it's data tag and data to an ofstream, and
 * can append data to the end of the variable once it has been written
 * to the file.
 */ 
class MatlabArray : public MatlabElement {
  friend class MatlabDataWriter; 
private:

  // Not allowed
  MatlabArray(const MatlabArray &rhs)
  {}
  const MatlabArray &operator=(const MatlabArray &rhs)
  {
    return *this;
  }

protected:
  MatlabElement me_flags_;
  MatlabElement me_dim_lens_;
  MatlabElement me_name_;
  MatlabElement me_data_;

  array_flags_t flags_; /**< Flags for this array. */ 
  string name_; /**< The name of this array */
  bool time_dim_; /**< True if this thing has a time dimension. */ 

  unsigned int num_dims_; /**< Number of dimensions */ 
  int32_t *dim_lengths_; /**< Dimension lengths */ 

  /**
   * Translate the given MATLAB data type into the appropriate MATLAB
   * array class.
   *
   * @param type the MATLAB data type
   */ 
  MATLAB_array_type get_array_class(MATLAB_data_type type);
  
public:
  MatlabArray(const char *name, 
              const vector<int> &dim_lens, 
              bool time_dim,
              MPI_Datatype type,
              bool complex); 

  virtual ~MatlabArray();

  /**
   * Append to the buffer. 
   *
   * @param num_bytes The number of bytes to append to the end of
   * this element. 
   * @param ptr The memory to get the data from. 
   */
  virtual void append_buffer(unsigned int num_bytes, const void *ptr);

  /**
   * Write the buffered data to disk
   *
   * @param stream the stream to write the data to
   */ 
  virtual void write_buffer(ostream &stream);

  /**
   * How many bytes (data and tag) will this variable write? This
   * must ask all of the children how big they are, and report the
   * sum plus the header. 
   */
  virtual unsigned int get_num_bytes();
};

/**
 * Writes data to a Matlab 5 binary data file. Data is aligned to 64
 * bit boundaries in the file. For matrix types, the number of bytes
 * in the tag includes these any padding bytes that may be required
 * for alignment; for other types the padding bytes are NOT included.
 *
 * \bug This only implements matrix output... it shouldn't be too much
 * of a stretch to implement everything, but we don't really need it.
 *
 * \bug All data is buffered until deinit() is called...
 */
class MatlabDataWriter : public DataWriter {
  
private:
protected:
  string var_name_; /**< For future comparison. */

  ofstream file_;
  
  unsigned int dim_len_;

  struct {
    char text[124];
    uint16_t version;
    uint8_t endian[2]; 
  } header_;

  map<string, MatlabArray *> vars_;

  bool buffered_; /**< Data is buffered in memory and only written out
                     at the end, when deinit() is called. */ 

  /**
   * Called by gather_data(), calls the above function. 
   */
  unsigned int write_data(unsigned int time_step, 
                          Variable &var, MPI_Datatype t, void *ptr, 
                          unsigned int len);

  /**
   * Called by gather_data() to indicate that all of the data
   * availble has been written. 
   */
  void end_data();

  /**
   * Set up the header for the file
   */ 
  void header_setup();

  /**
   * Add a numeric array with a given MATLAB datatype. 
   *
   * @param name The name of the variable (in the MATLAB workspace)
   * @param num_dimensions The number of dimensions of the array
   * @param type The MATLAB datatype of the array
   * @param complex True if the array contains elements which are
   * complex. 
   */
  void add_array(const char *name, 
                 const vector<int> &dimensions, 
                 MATLAB_data_type type,
                 bool complex);

public:
  MatlabDataWriter(int rank, int size);
  MatlabDataWriter(int rank, int size, const char *fn, Result &result);
  virtual ~MatlabDataWriter();

  /**
   * Initialize this object; open the file
   */
  virtual void init(const Grid &grid);

  /**
   * Deinit; close the file.
   */
  virtual void deinit(const Grid &grid);

  /**
   * Add a variable that we should know about. 
   * @param result describes the variable
   */
  virtual void add_variable(Result &result);

  /**
   * Write a test file
   */ 
  void test();

};

#endif // ACSII_DATA_WRITER_H
