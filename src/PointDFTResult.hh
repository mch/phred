#ifndef POINT_DFT_RESULT_H
#define POINT_DFT_RESULT_H

#include "Result.hh"
#include "TimeExcitation.hh"

/**
 * Outputs the DFT of the a field component at a point in space. 
 */
class PointDFTResult : public Result
{
private:
protected:
  field_t freq_start_; /**< First frequency in the range */
  field_t freq_stop_; /**< Last frequency in the range */
  unsigned int num_freqs_; /**< Number of frequencies in range */
  field_t freq_space_;

  field_t *result_; /**< Storage for the result. Interlevaved data;
                       freq, Real DFT value, Imag DFT value, etc */

  point_t point_; /**< The point to aquire the data at. */

public:
  PointDFTResult(field_t freq_start, field_t freq_stop, 
                 unsigned int num_freqs);
  ~PointDFTResult();

  /**
   * Set the point in global coordinates
   *
   * @param point
   */
  inline void set_point(point_t p)
  {
    point_ = p;
  }

  /**
   * Get the point in global coordinates
   *
   * @return point_t
   */
  inline point_t get_point()
  {
    return point_;
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
  virtual Data &get_result(Grid &grid, unsigned int time_step);

};

#endif // POINT_DFT_RESULT_H
