#ifndef DATA_WRITER_H
#define DATA_WRITER_H

#include "Types.hh"

#include <map>
#include <vector>

using namespace std;

class DataWriter;

/**
 * A variable which defines the name, dimensionality, and size of
 * the data to be output. 
 */
class Variable
{
private:
  Variable();

protected:
  string name_;

  DataWriter dw_; /**< Output goes here */

public:
  Variable(DataWriter &dw, string name)
    : name_(name), dw_(dw)
  {}

  virtual Variable() = 0;
};

/**
 * A variable for 4d (x,y,z,t) data
 */
class Variable_4d : public Variable
{
private:
  Variable_4d();
protected:
  region r_;

public:
  Variable_4d(DataWriter &dw, string name, region r)
    : Variable(dw, name), r_(r)
  {}
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
   * Add a 4d variable; x, y, z, and time! The x, y, and z axis are
   * bounded by region, but time is assumed to be infinitly long,
   * unless that's not supported. 
   *
   * @param name a variable name
   * @param r describes the region to output
   */
  Variable_4d* add_4d_variable(string name, region r);
};

#endif // DATA_WRITER_H
