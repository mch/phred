#ifndef EXCITE_TIME_H
#define EXCITE_TIME_H

#include "Excitation.hh"

/**
 * A type of excitation that depends only on the time step at which it
 * is being used. This class has a loop in it which calls the excite()
 * function of subclasses once and then applies it to every point the
 * region where it is to be applied.
 */
class TimeExcitation : public Excitation
{
protected:

public:
  TimeExcitation() {}
  virtual ~TimeExcitation() {}

  /**
   * This function is called to apply this excitation to the grid. 
   *
   * @param the grid to operate on
   * @param the time step we are on
   */
  virtual void excite(Grid &grid, unsigned int time_step);

  /**
   * This function is defined in subclasses and produces the source
   * value to be applied to every point which is excited.
   *
   * @param the grid to which the function is being applied. Only for
   * things such as delta_t, the result of this function is applied in
   * the excite() function (in order to factor out the loop). 
   *
   * @param the grid to which the function is being applied. Only for
   * things such as delta_t, the result of this function is applied in
   * the excite() function (in order to factor out the loop). 
   *
   * @param the time at which to apply the excitation
   */
  virtual field_t source_function(Grid &grid, unsigned int time_step) = 0;

};

#endif // EXCITE_TIME_H
