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

#ifndef DATA_H
#define DATA_H

#include <string>
#include <map>

#include <mpi.h>

using namespace std;

/**
 * A wrapper for a MPI derived data type, a pointer, and an int
 * indicating the number of items. This breaks OO principles slightly,
 * since we are exposing pointers to member data, but oh well. 
 *
 * Data object can contain multiple pointers, but they must all use
 * the same MPI datatype. An example use is for the BlockResult,
 * which can optionally return more than one field component. 
 */
class Data {
private:
protected:
  MPI_Datatype type_; /**< LOCAL datatype. For interperting data sent
                       * from a single node. */
  MPI_Datatype global_type_; /**< GLOBAL datatype. For interperting
                              * the data once it has been collected to
                              * the output node. */ 

  map<unsigned int, void *> ptrs_;
  unsigned int num_; /**< Number of items of type_ that can be
                        accessed by ptr_ */ 

public:
  Data()
    : type_(GRID_MPI_TYPE), global_type_(MPI_DATATYPE_NULL), num_(0)
  {
    //MPI_Type_contiguous(1, GRID_MPI_TYPE, &type_);
    //MPI_Type_commit(&type_);
  }

  ~Data()
  {}

  /**
   * Used to set the data type.
   * @param type MPI derived data type
   */
  inline void set_datatype(MPI_Datatype &type)
  {
    type_ = type;
  }

  /**
   * Get the data type
   */
  inline MPI_Datatype get_datatype() const
  {
    return type_;
  }

  /**
   * Used to set the global data type.
   *
   * @param type MPI derived data type
   */
  inline void set_global_datatype(MPI_Datatype &type)
  {
    global_type_ = type;
  }

  /**
   * Get the global data type. DataWriters should use this data type
   * to interpert the collected chuck of data that is to be
   * written. For example, one process may return a column of ints,
   * and another may return a column of floats, as reflected by each
   * node's local data type for this variable. This datatype specifies
   * how those two columns go together. 
   *
   * The more common case is a plane or block result, where this
   * datatype specifies the entire (global) region of data, but the
   * local data type species only the data chuck local to that node.
   */
  inline MPI_Datatype get_global_datatype() const
  {
    if (MPI_DATATYPE_NULL == global_type_)
      return type_;
    else
      return global_type_;
  }

  /**
   * Returns the number of pointers contained in this data
   * block. 
   */
  inline unsigned int get_num_ptrs() const 
  {
    return ptrs_.size();
  }

  /** 
   * Returns a pointer to some data. Data makes no guarantee about the
   * state of the pointer, because it is set by some other object. 
   *
   * @param ptr_num the pointer number to return, defaults to
   * zero. Must be less than the value returned by get_num_ptrs(). 
   */
  inline void *get_ptr(unsigned int ptr_num = 0) const
  {
    map<unsigned int, void *>::const_iterator iter = ptrs_.find(ptr_num);
    
    if (iter != ptrs_.end())
      return (*iter).second;
    else
      return 0;
  }

  /**
   * Returns the number of items (MPI Datatype data blocks) contained
   * in each pointer. 
   */
  inline unsigned int get_num() const
  {
    return num_;
  }

  /**
   * Sets a pointer. The optional parameter ptr_num determines which
   * pointer number to set. 
   *
   * @param ptr_num optional parameter the specifies the pointer
   * number to set. Defaults to zero. 
   */
  inline void set_ptr(void *ptr, unsigned int ptr_num = 0)
  {
    ptrs_[ptr_num] = ptr;
  }

  /**
   * Set the number of items. Setting to zero indicates no result. 
   */
  inline void set_num(unsigned int num)
  {
    num_ = num;
  }
};


/**
 * This is an extension of the Data class which allows for two data
 * pointers; a recieve (rx) and a transmit (tx). This is intended for
 * use when there is overlapping regions at subdomain boundaries and
 * data has to be transfered between ranks. 
 *
 * This is really just a convienence class. 
 */
class RxTxData : public Data
{
private:
  // Hide these
  inline void *get_ptr(unsigned int ptr_num) const
  { return Data::get_ptr(ptr_num); }

  inline void set_ptr(void *ptr, unsigned int ptr_num)
  { Data::set_ptr(ptr, ptr_num); }

  inline unsigned int get_num_ptrs() const
  { return 0; }
  
protected:
  FieldType field_type_; /**< The field type of the data being exchanged */

public:
  RxTxData()
    : field_type_(E)
  {}

  /**
   * Set the rx pointer
   */
  inline void set_rx_ptr(void *ptr)
  {
    set_ptr(ptr, 0);
  }

  /**
   * Get the rx pointer
   */
  inline void *get_rx_ptr() const
  {
    return get_ptr(0);
  }

  /**
   * Set the tx pointer
   */
  inline void set_tx_ptr(void *ptr)
  {
    set_ptr(ptr, 1);
  }

  /**
   * Get the tx pointer
   */
  inline void *get_tx_ptr() const
  {
    return get_ptr(1);
  }

  /**
   * Set the field type (E or H)
   */
  inline void set_field_type(FieldType t)
  {
    field_type_ = t;
  }

  /**
   * Returns the field type of this data
   */
  inline FieldType get_field_type()
  {
    return field_type_;
  }
  
};

#endif // DATA_H
