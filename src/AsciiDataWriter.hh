#ifndef ACSII_DATA_WRITER_H
#define ACSII_DATA_WRITER_H

#include "DataWriter.hh"
#include <fstream>
#include <exception>

using namespace std;

/**
 * Writes data to a text file that can be loaded by MATLAB. 
 */
class AsciiDataWriter : public DataWriter {
  
private:
protected:
  string filename_;
  string var_name_; /**< For future comparison. */

  ofstream file_;
  
  unsigned int dim_len_;

  /**
   * Does recursive writing of packed data. Increments the pointer
   * after writing each value, and returns it. 
   */
  void *write_data(MPI_Datatype t, void *ptr);

public:
  AsciiDataWriter(int rank, int size);
  ~AsciiDataWriter();

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

  /**
   * Add a variable that we should know about. Throws an exception if
   * you try to add more than one. This call currently only supports
   * one dimensional data as well. 
   *
   * @param result describes the variable
   */
  void add_variable(Result &result);

  /**
   * Write the Data produced by the Result object to a file. 
   *
   * @param data a Data object containing the data to handle
   */
  void handle_data(unsigned int time_step, Data &data);

};

#endif // ACSII_DATA_WRITER_H
