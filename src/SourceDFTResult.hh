#ifndef SOURCE_DFT_RESULT_H
#define SOURCE_DFT_RESULT_H

#include "Result.hh"
#include "SourceFunction.hh"

/**
 * Outputs the DFT of a source function. Only really applies to
 * excitations that are applied at a single point in space.
 */
class SourceDFTResult : public Result
{
private:
protected:
  SourceFunction &te_; /**< The dft excitation to save. */

  field_t freq_start_; /**< First frequency in the range */
  field_t freq_stop_; /**< Last frequency in the range */
  unsigned int num_freqs_; /**< Number of frequencies in range */
  field_t freq_space_;

  field_t *result_; /**< Storage for the result. Interlevaved data;
                       freq, Real DFT value, Imag DFT value, etc */

public:
  SourceDFTResult(SourceFunction &te);
  SourceDFTResult(SourceFunction &te, field_t freq_start,
                  field_t freq_stop, unsigned int num_freqs);
  ~SourceDFTResult();

  /**
   * Set the source or TimeExcitation to use.
   * @param a reference to a TimeExcitation object
   */ 
  void set_excitation(const SourceFunction &te);

  /**
   * Set the start frequency of the range
   */
  inline void set_freq_start(field_t fs)
  {
    freq_start_ = fs;
  }

  /**
   * Set the end frequency of the range
   */
  inline void set_freq_stop(field_t fs)
  {
    freq_stop_ = fs;
  }

  /**
   * Set the number of frequencies of in the range
   */
  inline void set_num_freq(unsigned int nf)
  {
    num_freqs_ = nf;
  }

  /**
   * Get the start frequency of the range
   */
  inline field_t get_freq_start()
  {
    return freq_start_;
  }

  /**
   * Get the end frequency of the range
   */
  inline field_t get_freq_stop()
  {
    return freq_stop_;
  }

  /**
   * Get the number of frequencies of in the range
   */
  inline unsigned int get_num_freq()
  {
    return num_freqs_;
  }

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
  virtual Data &get_result(const Grid &grid, unsigned int time_step);

  /**
   * Setup the result, allocate memory, etc. Called just before the
   * simulation starts. 
   */
  void init(const Grid &grid);
  
  /**
   * Deallocates memory. 
   */
  void deinit(const Grid &grid);

};

#endif // SOURCE_DFT_RESULT_H
