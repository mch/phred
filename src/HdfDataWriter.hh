#ifndef HDF_DATA_WRITER_H
#define HDF_DATA_WRITER_H

#include "DataWriter.hh"
#include "config.h"
#include "Exceptions.hh"

#ifdef USE_HDF4

#include <map>
#include <vector>

using namespace std;

// I don't want the class to disappear when using Python scripts, I
// just want it to throw an exception when you try to use it and you
// can't.  Then you can look for the exception in your script and try
// a different data writer.

/**
 * HDF4 data writer; collects all data to be written to one rank and
 * saves it to disk.
 *
 * \bug IMPLEMENT ME!
 */
class HdfDataWriter : public DataWriter {
private:
protected:
  string filename_;


public:
  HdfDataWriter(int rank, int size)
    : DataWriter(rank, size)
  {
    throw NoHDF4Exception(); 
  }

  ~HdfDataWriter()
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

class HdfDataWriter : public DataWriter {
public:
  HdfDataWriter(int rank, int size)
    : DataWriter(rank, size)
  {
    throw NoHDF4Exception(); 
  }

  ~HdfDataWriter()
  {}

  void init()
  {}

  void deinit()
  {}

  void add_variable(Result &result)
  {}

};

#endif // USE_HDF4

#endif // HDF_DATA_WRITER_H
