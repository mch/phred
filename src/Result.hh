#ifndef RESULT_H
#define RESULT_H

#include <string>
#include "Types.hh"
#include "Grid.hh"

using namespace std;

/**
 * A wrapper for a MPI derived data type, a pointer, and an int
 * indicating the number of items. 
 */
class Data {
private:
protected:
  MPI_Datatype type_;
  void *ptr_;
  unsigned int num_;

public:
  Data();
  ~Data();

  /**
   * Used to set the data type.
   */
  inline void set_datatype(MPI_Datatype &type)
  {
    type_ = type;
  }

  /**
   * Get the data type
   */
  inline MPI_Datatype get_datatype()
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
 */
class Result
{
private:
protected:
  string var_name_; /**< Variable name */

public:
  Result() {}
  virtual ~Result() = 0;

  /**
   * Looks at the grid and produces output
   *
   * @param grid a reference to a Grid object
   * @return a data object, which contains an MPI derived data type,
   * a pointer, and the number of items in the result. 
   */
  virtual Data get_result(Grid &grid) = 0;

};

#endif // RESULT_H
