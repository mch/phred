/* 
   Phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004 Matt Hughes <mhughe@uvic.ca>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#ifndef DFT_RESULT_H
#define DFT_RESULT_H

#include "Result.hh"
#include "../Interval.hh"

/**
 * This is a little helper class that has all the common stuff that is
 * required to do DFT's in results. 
 */ 
class DFTResult : public Result
{
public:
  DFTResult();

  /**
   * Constructor. Sets the frequency parameters for this result.
   *
   * @param freq_start first frequency in the range
   * @param freq_stop last frequency in the range
   * @param num_freqs total number of frequencies to report data for
   */
  DFTResult(field_t freq_start, field_t freq_stop, 
            unsigned int num_freqs);

  virtual ~DFTResult();

  /**
   * Set the frequency parameters for this result.
   *
   * @param freq_start first frequency in the range
   * @param freq_stop last frequency in the range
   * @param num_freqs total number of frequencies to report data for
   */ 
  void set_freq(field_t freq_start, field_t freq_stop, 
                unsigned int num_freqs)
  { frequencies_.set_params(freq_start, freq_stop, num_freqs); }

  /**
   * Get the start frequency of the range
   */
  inline field_t get_freq_start() const
  {
    return frequencies_.get_start();
  }

  /**
   * Get the end frequency of the range
   */
  inline field_t get_freq_stop() const
  {
    return frequencies_.get_end();
  }

  /**
   * Get the number of frequencies of in the range
   */
  inline unsigned int get_num_freq() const
  {
    return frequencies_.length();
  }

  /**
   * Print a string representation to an ostream.
   */
  virtual ostream& to_string(ostream &os) const;

protected:

  Interval<field_t> frequencies_;

private:
};

#endif
