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
  MPI_Datatype type_;
  map<unsigned int, void *> ptrs_;
  unsigned int num_; /**< Number of items of type_ that can be
                        accessed by ptr_ */ 
  string var_name_; /**< The name of the variable this data belongs
                       to. */ 

public:
  Data(string var_name) 
  {
    // A single field_t by default. 
    MPI_Type_contiguous(1, GRID_MPI_TYPE, &type_);
    MPI_Type_commit(&type_);
    var_name_ = var_name;
  }

  Data()
  {
    MPI_Type_contiguous(1, GRID_MPI_TYPE, &type_);
    MPI_Type_commit(&type_);
  }

  ~Data()
  {}

  /**
   * Get the variable name for this data. 
   */
  inline string &get_var_name()
  {
    return var_name_;
  }

  /**
   * Set the variable name for this data.
   */
  inline void set_var_name(const string &var_name)
  {
    var_name_ = var_name;
  }

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
  inline MPI_Datatype &get_datatype()
  {
    return type_;
  }

  /**
   * Returns the number of pointers contained in this data
   * block. 
   */
  inline unsigned int get_num_ptrs()
  {
    return ptrs_.size();
  }

  /** 
   * Returns a pointer to some data. Data makes no guarantee about the
   * state of the point, because it is set by some other object. 
   *
   * @param ptr_num the pointer number to return, defaults to
   * zero. Must be less than the value returned by get_num_ptrs(). 
   */
  inline void *get_ptr(unsigned int ptr_num = 0)
  {
    map<unsigned int, void *>::iterator iter = ptrs_.find(ptr_num);
    
    if (iter != ptrs_.end())
      return (*iter).second;
    else
      return 0;
  }

  /**
   * Returns the number of items (MPI Datatype data blocks) contained
   * in each pointer. 
   */
  inline unsigned int get_num()
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
  inline void *get_ptr(unsigned int ptr_num)
  { return Data::get_ptr(ptr_num); }

  inline void set_ptr(void *ptr, unsigned int ptr_num)
  { Data::set_ptr(ptr, ptr_num); }

  inline unsigned int get_num_ptrs()
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
  inline void *get_rx_ptr()
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
  inline void *get_tx_ptr()
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
