#ifndef DATA_WRITER_H
#define DATA_WRITER_H

#include "Types.hh"

#include <map>
#include <vector>

using namespace std;

/**
 * An abstract base class that can be subclassed to create objects
 * which save data to disk, or perform other output
 * functions. DataWriter's have a list of variables, although some
 * writers may only be able to support one. Variables have dimensions,
 * units, and names. DataWriter's must support at least one
 * dimensional data.
 *
 * \bug add_variable doesn't check if a varable by that name already exists!
 */
class DataWriter 
{
private:
public:

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

    virtual ~Variable() = 0;
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

    /** 
     * Implementing classes must implment this to do the actual saving. 
     */
    virtual void write_point(unsigned int x, unsigned int y,
                             unsigned int z, field_t val) = 0;

  public:
    Variable_4d(DataWriter &dw, string name, region r)
      : Variable(dw, name), r_(r)
    {}

    virtual ~Variable_4d() = 0;

    /**
     * Save a point of data. Simple. Slow. Sad. If this DataWriter is
     * not on rank 0, the data is transmitted to rank 0. And then the
     * write_point function is called by rank 0. 
     *
     * @param x x coord
     * @param y y coord
     * @param z z coord
     * @param val value
     */ 
    void save_point(unsigned int x, unsigned int y,
                    unsigned int z, field_t val);
  };

protected:
  unsigned int max_vars_; /**< Maximum number of variable for this
                             file */

  map<string, *Variable> variables_;
  
  int rank_;  /**< MPI Rank of this process */
  int size_; /**< Number of processes in MPI comm. */

private:
  DataWriter()
  {}

public:
  DataWriter(int rank, int size) 
    : max_vars_(0), rank_(rank), size_(size) 
  {}

  virtual ~DataWriter() = 0;

  /**
   * Initialize this object
   */
  virtual void init() = 0;

  /**
   * Deinit
   */
  virtual void deinit() = 0;

  /**
   * Returns the rank
   */
  inline int get_rank()
  {
    return rank_;
  }

  /**
   * Returns the size
   *
   * @return an int; number of MPI processes
   */
  inline int get_size()
  {
    return size_;
  }

  /**
   * Add a 4d variable; x, y, z, and time! The x, y, and z axis are
   * bounded by region, but time is assumed to be infinitly long,
   * unless that's not supported. 
   *
   * @param name a variable name
   * @param r describes the region to output
   */
  Variable_4d& add_4d_variable(string name, region r);

};

#endif // DATA_WRITER_H
