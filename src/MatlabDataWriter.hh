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
}

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
  uint8_t undefined[2];
  union flags {
    uint8_t byte;
    struct bits {
      unsigned bits : 8;
    }
  }
  uint8_t array_class;
  uint8_t undefined[4];
} array_flags_t;

/**
 * Represents a data element in a MATLAB file. 
 */
class MatlabElement {
  friend class MatlabMatrix;
  friend class MatlabDataWriter;
private:
  MatlabElement();
  MatlabElement(const MatlabElement &rhs);
  const MatlabElement &operator=(const MatlabElement &rhs);

protected:
  data_tag_t tag_; /**< Type and size of data */ 
  unsigned char num_pad_bytes_; /**< Number of pad bytes on the end of
                                   the data to make it align to 64 bit
                                   boundaries. */  
  
  unsigned int file_offset_; /**< Position (in bytes) in the file of
                                the start of the tag for this element */ 
public:

  /**
   * How many bytes (data and tag) will this variable write? 
   */
  virtual unsigned int get_num_bytes();

  /**
   * Set the MATLAB data type for this element. 
   *
   * @param type The MATLAB datatype
   */
  inline void set_type(MATLAB_data_type type)
  {
    tag_.datatype = type;
    tag_.num_bytes = 0;
  }

  /**
   * Writes some data to the given ostream
   *
   * @param stream The output stream to write to. 
   * @param num_bytes The number of bytes to append to the end of
   * this element. 
   * @param ptr The memory to get the data from. 
   */
  virtual void write_data(ostream &stream, unsigned int num_bytes,
                          const void *ptr);
  
  /**
   * Updates the file offset so that we know where to find the tag in
   * the file. 
   *
   * @param bytes The number of bytes to add to the current offset. 
   */ 
  virtual void update_file_offset(int bytes);
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
  MatlabArray(); 
  MatlabArray(const MatlabVariable &rhs);
  const MatlabArray &operator=(const MatlabVariable &rhs);

protected:
  vector <MatlabElement *> children_;

public:
  /**
   * Construct a numeric array with a given MATLAB datatype. 
   *
   * @param name The name of the variable (in the MATLAB workspace)
   * @param num_dimensions The number of dimensions of the array
   * @param type The MATLAB datatype of the array
   * @param complex True if the array contains elements which are
   * complex. 
   */
  static MatlabArray build_array(const char *name, 
                                 unsigned int num_dimensions, 
                                 MATLAB_data_type type,
                                 bool complex);

  //static MatlabArray build_char_array(const char *name);

  /**
   * Writes this array to the given ostream
   *
   * @param stream The output stream to write to. 
   * @param num_bytes The number of bytes to append to the end of
   * this element. 
   * @param ptr The memory to get the data from. 
   */
  void write_data(ostream &stream, unsigned int num_bytes,
                          const void *ptr);
  
  /**
   * How many bytes (data and tag) will this variable write? This
   * must ask all of the children how big they are, and report the
   * sum plus the header. 
   */
  unsigned int get_num_bytes();
};

/**
 * Writes data to a Matlab 5 binary data file. Data is aligned to 64
 * bit boundaries in the file. For matrix types, the number of bytes
 * in the tag includes these any padding bytes that may be required
 * for alignment; for other types the padding bytes are NOT included. 
 */
class MatlabDataWriter : public DataWriter {
  
private:
protected:
  string var_name_; /**< For future comparison. */

  ofstream file_;
  
  unsigned int dim_len_;

  struct header_ {
    char text[124];
    uint16_t version;
    union endian {
      uint16_t endian_16; 
      char bytes[2];
    }
  }

  map<string, MatlabVariable *> vars_;

  /**
   * Write an array of data
   * @param ptr location of start of array
   * @param len number of items in array
   * @return the advanced ptr
   */
  template<class T>
  unsigned int write_array(void *ptr, unsigned int len);

  /**
   * Does recursive writing of packed data. This function will only be
   * called on rank 0.
   */
  unsigned int write_data(MPI_Datatype t, void *ptr, 
                          unsigned int len);

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

public:
  MatlabDataWriter(int rank, int size);
  MatlabDataWriter(int rank, int size, const char *fn, Result &result);
  ~MatlabDataWriter();

  /**
   * Initialize this object; open the file
   */
  void init();

  /**
   * Deinit; close the file.
   */
  void deinit();

  /**
   * Add a variable that we should know about. 
   * @param result describes the variable
   */
  void add_variable(Result &result);

};

#endif // ACSII_DATA_WRITER_H
