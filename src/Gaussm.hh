/** \class Gaussm
 * \brief A modulated gaussian excitation
 *
 * A Gaussian function modulated by a sine wave. 
 */

#ifndef EXCITE_GAUSSM_H
#define EXCITE_GAUSSM_H

#include "SourceFunction.hh"
#include "Constants.hh"

#include <math.h>

class Gaussm : public SourceFunction
{
protected:
  field_t alpha_;
  field_t deltaf_;
  field_t f0_;
  
public:
  Gaussm();
  virtual ~Gaussm();

  /**
   * Set the parameters of the modulated Gauss function.
   *
   * @param alpha scale factor (defaults to 1)
   * @param deltaf the frequency width of the pulse (defaults to 50 MHz)
   * @param f0 the centre frequency of the pulse (defaults to 1 GHz)
   */
  void set_parameters(field_t alpha, field_t deltaf, field_t f0);

  /**
   * Get the current alpha value
   *
   * @return alpha
   */
  field_t get_alpha();

  /**
   * Get the current frequency range
   *
   * @return frequency range
   */
  field_t get_deltaf();

  /**
   * Get the current centre frequency
   *
   * @return centre frequency
   */
  field_t get_f0();

  /**
   * Produces a modulated gauss function. 
   *
   * @param grid the grid to which the function is being applied. Only for
   * things such as delta_t, the result of this function is applied in
   * the excite() function (in order to factor out the loop). 
   *
   * @param time_step the time at which to apply the excitation
   *
   * @return the value of the excitation
   */
  field_t source_function(const Grid &grid, unsigned int time_step);

};

#endif // EXCITE_GAUSSM_H
