// A modulated gaussian excitation

#ifndef EXCITE_GAUSSM_H
#define EXCITE_GAUSSM_H

#include "TimeExcitation.hh"
#include "Constants.hh"

#include <math.h>

class Gaussm : public TimeExcitation
{
protected:
  field_t alpha_;
  field_t deltaf_;
  field_t f0_;
  
public:
  Gaussm();
  ~Gaussm();

  /**
   * Set the parameters of the modulated Gauss function.
   *
   * @param scale factor (defaults to 1)
   * @param the frequency width of the pulse (defaults to 50 MHz)
   * @param the centre frequency of the pulse (defaults to 1 GHz)
   */
  void set_parameters(field_t alpha, field_t deltaf, field_t f0);

  /**
   * Get the current alpha value
   */
  field_t get_alpha();

  /**
   * Get the current frequency range
   */
  field_t get_deltaf();

  /**
   * Get the current centre frequency
   */
  field_t get_f0();

  /**
   * Produces a modulated gauss function. 
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
  field_t source_function(Grid &grid, unsigned int time_step);

};

#endif // EXCITE_GAUSSM_H
