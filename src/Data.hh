#ifndef DATA_H
#define DATA_H

#include <string>
#include <mpi.h>

using namespace std;

/**
 * A wrapper for a MPI derived data type, a pointer, and an int
 * indicating the number of items. This breaks OO principles slightly,
 * since we are exposing pointers to member data, but oh well. 
 */
class Data {
private:
protected:
  MPI_Datatype type_;
  void *ptr_;
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
   * Get the pointer
   */
  inline void *get_ptr()
  {
    return ptr_;
  }

  /**
   * get the number of items
   */
  inline unsigned int get_num()
  {
    return num_;
  }

  /**
   * Set the pointer.
   */
  inline void set_ptr(void *ptr)
  {
    ptr_ = ptr;
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
 */
class RxTxData : public Data
{
private:
  // Hide these
  void *get_ptr()
  { return 0; }

  void set_ptr(void *)
  {}
  
protected:
  void *rx_ptr_;
  void *tx_ptr_;
public:
  /**
   * Set the rx pointer
   */
  inline void set_rx_ptr(void *ptr)
  {
    rx_ptr_ = ptr;
  }

  /**
   * Get the rx pointer
   */
  inline void *get_rx_ptr()
  {
    return rx_ptr_;
  }

  /**
   * Set the tx pointer
   */
  inline void set_tx_ptr(void *ptr)
  {
    tx_ptr_ = ptr;
  }

  /**
   * Get the tx pointer
   */
  inline void *get_tx_ptr()
  {
    return tx_ptr_;
  }

};

#endif // DATA_H
