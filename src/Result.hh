#ifndef RESULT_H
#define RESULT_H

#include <string>
#include "Types.hh"
#include "Grid.hh"

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
  unsigned int num_bytes_; /**< Number of bytes per item */
  string var_name_; /**< The name of the variable this data belongs
                       to. */ 

public:
  Data(string var_name) 
  {
    // A single field_t by default. 
    MPI_Type_contiguous(1, GRID_MPI_TYPE, &type_);
    MPI_Type_commit(&type_);
    var_name_ = var_name;
    num_bytes_ = sizeof(field_t);
  }

  Data()
  {
    MPI_Type_contiguous(1, GRID_MPI_TYPE, &type_);
    MPI_Type_commit(&type_);
    num_bytes_ = sizeof(field_t);
  }

  ~Data()
  {}

  /**
   * Get the number of bytes per item
   */
  inline unsigned int get_num_bytes()
  {
    return num_bytes_;
  }

  /**
   * Set the number of bytes per item.
   */
  inline void set_num_bytes(unsigned int nb)
  {
    num_bytes_ = nb;
  }

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
   * Set the number of items
   */
  inline void set_num(unsigned int num)
  {
    num_ = num;
  }
};

/**
 * An abstract base class used to implement results. This class
 * recieves a reference to a grid and a reference to a DataWriter. It
 * takes data from the grid, performs calculations, and then calls
 * functions in the DataWriter to save the data. 
 *
 * The method get_result() MUST be implemented. Is intended to be called
 * after every the updates for each time step have been computed. 
 *
 * get_result(), returns a structure containing an MPI derived
 * datatype describing the data, a pointer to the data, and the number
 * of elements in the result. If there is no data to be written, then
 * the number of elements is zero. 
 *
 * The dimensions defined by the subclasses is NOT to include
 * time. Since this request gets called on once every time step, time
 * is implied. 
 */
class Result
{
private:
  Result(const Result &rhs);
  const Result &operator=(const Result &rhs);

protected:
  string var_name_; /**< Variable name */
  Data data_;

  unsigned int num_dims_; /**< Number of dimensions */
  unsigned int *dim_lens_; /**< Dimension lengths */

  string dw_name_;
  //DataWriter *dw_;
public:
  Result() 
    : num_dims_(0), dim_lens_(0) //, dw_(0)
  {}

  virtual ~Result()
  {}

  /**
   * Set the data writer name (intended where we read that from a
   * config file and set the actual data writer later). 
   */
  inline void set_dw_name(const string &dw)
  {
    dw_name_ = dw;
  }

  /**
   * Get the data writer name
   */
  inline const string &get_dw_name()
  {
    return dw_name_;
  }

  /**
   * Looks at the grid and produces output
   *
   * @param grid a reference to a Grid object
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  virtual Data &get_result(Grid &grid, unsigned int time_step) = 0;

  /**
   * Set the name of the variable. 
   * @param name a string with the name
   */
  inline void set_name(const string &name)
  {
    var_name_ = name;
    data_.set_var_name(name);
  }

  /**
   * Return the name of the variable. 
   * @return a string with the name in it. 
   */
  inline const string &get_name()
  {
    return var_name_;
  }

  /**
   * Called to perform any initialization that may be required. Does
   * nothing by default. 
   */
  virtual void init()
  {}

  /** 
   * Returns the number of dimensions
   * @return unsigned int, num dims
   */
  inline unsigned int get_num_dims()
  {
    return num_dims_;
  }

  /**
   * Returns the lengths of the dimensions
   * @return a pointer to the lengths of the dimensions, with the 
   * length returned by get_num_dims()
   */
  inline unsigned int *get_dim_lengths()
  {
    return dim_lens_;
  }
};

#endif // RESULT_H
