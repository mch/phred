#ifndef EXCITE_TIME_H
#define EXCITE_TIME_H

#include "Excitation.hh"

/**
 * A type of excitation that depends on the time step at which it is
 * being used, and the point in space where it is being applied. This
 * class has a loop in it which calls the source_function() of
 * subclasses once and then applies it to every point the region where
 * it is to be applied.
 */
class SpaceExcitation : public Excitation
{
protected:

public:
  SpaceExcitation() {}
  virtual void ~SpaceExcitation() = 0;

  /**
   * This function is called to apply this excitation to the grid. 
   *
   * @param grid the grid to operate on
   * @param time_step the time step we are on
   * @param type the field type we are allowed to excite. 
   */
  virtual void excite(Grid &grid, unsigned int time_step,
                      FieldType type);

  /**
   * This function is defined in subclasses and produces the source
   * value for the given point in space and time. 
   *
   * @param grid the grid to which the function is being applied. Only for
   * things such as delta_t, the result of this function is applied in
   * the excite() function (in order to factor out the loop). 
   *
   * @param time_step the time at which to apply the excitation
   * @param x the x coordinate at which to apply the excitation. 
   * @param y the y coordinate at which to apply the excitation. 
   * @param z the z coordinate at which to apply the excitation. 
   */
  virtual field_t source_function(Grid &grid, unsigned int time_step, 
                                  unsigned int x, unsigned int y, 
                                  unsigned int z) = 0;

};

#endif // EXCITE_TIME_H
