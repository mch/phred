#ifndef SOURCE_DFT_RESULT_H
#define SOURCE_DFT_RESULT_H

#include "Result.hh"
#include "TimeExcitation.hh"

/**
 * Outputs the DFT of a source function at the end of it's run. 
 */
class SourceDFTResult : public Result
{
private:
protected:
  TimeExcitation &te_; /**< The dft excitation to save. */

  field_t freq_start_; /**< First frequency in the range */
  field_t freq_stop_; /**< Last frequency in the range */
  unsigned int num_freqs_; /**< Number of frequencies in range */
  field_t freq_space_;

  field_t *result_; /**< Storage for the result. Interlevaved data;
                       freq, Real DFT value, Imag DFT value, etc */

public:
  SourceDFTResult(TimeExcitation &te, field_t freq_start,
                  field_t freq_stop, unsigned int num_freqs);
  ~SourceDFTResult();

  /**
   * Set the source or TimeExcitation to use.
   * @param a reference to a TimeExcitation object
   */ 
  void set_excitation(const TimeExcitation &te);

  /**
   * Produces a running DFT from the source, considering only the time
   * step. Output is only available at time_end. 
   *
   * @param grid a reference to a Grid object
   * @param time_step current time step
   * @return a reference to a data object, which contains an MPI
   * derived data type, a pointer, and the number of items in the
   * result.
   */
  virtual Data &get_result(Grid &grid, unsigned int time_step);

};

#endif // SOURCE_DFT_RESULT_H
