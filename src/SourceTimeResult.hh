#ifndef SOURCE_TIME_RESULT_H
#define SOURCE_TIME_RESULT_H

#include "Result.hh"
#include "SourceFunction.hh"

/**
 * Outputs a source function at the timestep.
 */
class SourceTimeResult : public Result
{
private:
protected:
  SourceFunction &te_; /**< The time excitation to save. */

  field_t result_[2]; /**< Storage for the result. */

public:
  SourceTimeResult(SourceFunction &te);
  ~SourceTimeResult();

  /**
   * Set the source or SourceFunction to use.
   * @param a reference to a SourceFunction object
   */ 
  void set_excitation(const SourceFunction &te);

  /**
   * Produces output from the source, considering only the time step. 
   *
   * @param grid a reference to a Grid object
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  virtual Data &get_result(Grid &grid, unsigned int time_step);

};

#endif // SOURCE_TIME_RESULT_H
