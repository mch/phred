#ifndef ACSII_DATA_WRITER_H
#define ACSII_DATA_WRITER_H

#include "DataWriter.hh"

/**
 * Implements a

/**
 * Writes data to a text file that can be loaded by MATLAB. 
 */
class AsciiDataWriter : public DataWriter {
  
private:
protected:
  string filename_;

  bool open_;

public:
  class Variable_4d : public DataWriter::Variable_4d {
    void write_point(unsigned int x, unsigned int y,
                     unsigned int z, field_t val);
  };
  
  AsciiDataWriter(int rank, int size)
    : DataWriter(rank, size), open_(false)
  {}

  /**
   * Set the filename
   *
   * @param filename
   */
  inline void set_filename(string filename)
  {
    filename_ = filename;
  }

  /**
   * Returns the filename
   *
   * @return a string with the filename
   */
  inline string get_filename(string filename)
  {
    return filename_;
  }

  /**
   * Initialize this object; open the file
   */
  void init();

  /**
   * Deinit; close the file.
   */
  void deinit();

};

#endif // ACSII_DATA_WRITER_H
