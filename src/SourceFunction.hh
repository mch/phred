#ifndef SOURCE_FUNCTION_H
#define SOURCE_FUNCTION_H

#include "Types.hh"
#include "Grid.hh"

/**
 * This class defines a source_function which can be computed at some
 * timestep on some grid.
 *
 * Subclassess of this are intended to be template parameters to
 * TimeExcitation and SpaceTimeExcitation. 
 */
class SourceFunction
{
private:
protected:
public:
  SourceFunction() 
  {}
  virtual ~SourceFunction() 
  {}

  /**
   * This function is defined in subclasses and produces the source
   * value to be applied to every point which is excited.
   *
   * @param grid the grid to which the function is being applied. Only for
   * things such as delta_t, the result of this function is applied in
   * the excite() function (in order to factor out the loop). 
   *
   * @param time_step the time at which to apply the excitation
   */
  virtual field_t source_function(Grid &grid, unsigned int time_step) = 0;

};

#endif // SOURCE_FUNCTION_H
