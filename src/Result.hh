#ifndef RESULT_H
#define RESULT_H

#include <string>
#include "Types.hh"
#include "DataWriter.hh"
#include "Grid.hh"

using namespace std;

/**
 * An abstract base class used to implement results. This class
 * recieves a reference to a grid and a reference to a DataWriter. It
 * takes data from the grid, performs calculations, and then calls
 * functions in the DataWriter to save the data. 
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
   * @param writer a reference to a DataWriter object
   */
  virtual void produce_output(Grid &grid, DataWriter &dw) = 0;

};

#endif // RESULT_H
