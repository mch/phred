#ifndef HDF5_DATA_WRITER_H
#define HDF5_DATA_WRITER_H

#include "DataWriter.hh"
#include "config.h"
#include "Exceptions.hh"

#ifdef USE_HDF5

#include <map>
#include <vector>

using namespace std;

// I don't want the class to disappear when using Python scripts, I
// just want it to throw an exception when you try to use it and you
// can't.  Then you can look for the exception in your script and try
// a different data writer.

/**
 * HDF5 data writer; collects all data to be written to one rank and
 * saves it to disk.
 *
 * \bug IMPLEMENT ME!
 */
class Hdf5DataWriter : public DataWriter {
private:
protected:
  string filename_;


public:
  Hdf5DataWriter(int rank, int size)
    : DataWriter(rank, size)
  {
    throw NoHDF5Exception(); 
  }

  ~Hdf5DataWriter()
  {}

  /**
   * Set the filename to write to. 
   *
   * @param filename 
   */
  inline void set_filename(string filename)
  {
    filename_ = filename;
  }

  /**
   * Initialize this object; open the file. If the file is already
   * open, this simply returns. 
   */
  virtual void init();

  /**
   * Deinit; close the file. If the file is closed, this simply
   * returns. 
   */
  virtual void deinit();

  /**
   * Add a result that this data writer will have to know how to
   * handle. The file must be open (init() has to have been called)
   * before using this function. 
   *
   * @param result describes the result
   */
  virtual void add_variable(Result &result);

};

#else

class Hdf5DataWriter : public DataWriter {
public:
  Hdf5DataWriter(int rank, int size)
    : DataWriter(rank, size)
  {
    throw NoHDF5Exception(); 
  }

  ~Hdf5DataWriter()
  {}

  void init()
  {}

  void deinit()
  {}

  void add_variable(Result &result)
  {}

};

#endif // USE_HDF5

#endif // HDF5_DATA_WRITER_H
