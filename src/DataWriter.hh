#ifndef DATA_WRITER_H
#define DATA_WRITER_H

#include "Types.hh"

#include <map>
#include <vector>

using namespace std;

/** 
 * A Dimension
 */
class Dimension
{
protected:
  string label_;
  string units_;
  unsigned int len_;
};

/**
 * A variable which defines the name, dimensionality, and size of
 * the data to be output. 
 */
class Variable
{
protected:
  string name_;
  vector<*Dimension> dimensions_;

public:
  /**
   * Add dimension
   */
  Dimension *add_dimension();
};

/**
 * An abstract base class that can be subclassed to create objects
 * which save data to disk, or perform other output
 * functions. DataWriter's have a list of variables, although some
 * writers may only be able to support one. Variables have dimensions,
 * units, and names. DataWriter's must support at least one
 * dimensional data.
 */
template <class T>
class DataWriter 
{
private:
protected:
    
  map<string, *Variable> variables_;
  
public:

  DataWriter() {}
  virtual ~DataWriter() = 0;

  /**
   * Do something with a chunk of data.
   *
   * @param len number of items
   * @param data data to save
   */ 
  virtual void save_data(unsigned int len, T data) = 0;

  /**
   * Add a variable
   */
  Variable* add_variable(string name);
};

#endif // DATA_WRITER_H
