/* 
   Phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004-2005 Matt Hughes <mhughe@uvic.ca>

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


/** \class Excitation
 * \brief A base class for excitations. The default implementation
 * applies a source function to every point in a given
 * region. Subclasses may implement windowing functions or in other
 * ways depend on the position in the grid where the excitation is
 * being applied. 
 */

#ifndef EXCITE_H
#define EXCITE_H

#include "../config.h"

#include "../Types.hh"
#include "../Grid.hh"
#include "../Signals/Signal.hh"
#include "../LifeCycle.hh"

class Excitation : public LifeCycle
{
protected:
  /**
   * Region to apply the signal to (in the local grid)
   */
  shared_ptr<Block> region_;

  /**
   * Field polarization to excite, x, y, z
   */
  field_t polarization_[3];

  /**
   * Field to excite, E or H
   */
  FieldType type_;

  /**
   * True if this is a soft signal (defaults to false). Signal
   * signals are added to the field rather than the field being set
   * to them.
   */
  bool soft_;

  shared_ptr<Signal> sf_; /**< Signal function object to apply */

  shared_ptr<CSGBox> box_; /**< The box representing the area where
                              the excitation is to be applied. */
  
  float time_start_;
  float time_stop_;
  float time_offset_;

public:

  /**
   * Constructs a Excitation object which can apply a signal
   * function which depends only on time to a region of a
   * grid. Creates a copy of the signal function object to use, so you
   * don't have to hold onto a reference.
   *
   * @param sf Signal object to use as an excitation. 
   */
  Excitation(shared_ptr<Signal> sf);

  virtual ~Excitation();

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
   * Set the CSGBox inside which the excitation is to be applied. 
   *
   * @param box the CSGBox to apply the excitation to. 
   */
  virtual void set_region(shared_ptr<CSGBox> box);

  /**
   * Set the polarization vector
   *
   * @param x the x component
   * @param y the y component
   * @param z the z component
   */
  void set_polarization(field_t x, field_t y, field_t z);

  /**
   * Set the field to excite. Either E or H. 
   *
   * @param t FieldType
   */
  void set_type(FieldType t);

  /**
   * Returns the field type.
   */
  inline FieldType get_type()
  {
    return type_;
  }

  /**
   * Tells if this source is "soft" or not. Soft sources are
   * additive. 
   */
  inline bool get_soft()
  {
    return soft_;
  }

  /**
   * Set the softness of the source.
   * @param soft true if this source should be soft.
   */
  inline void set_soft(bool soft)
  {
    soft_ = soft;
  }

  /**
   * Verify that the excitation region is entierly inside the global grid. 
   */
  virtual void init(const Grid &grid);

  /**
   * Set the time window over which the excitation is applied, and the
   * time offset to be applied to the signal. All values are in real
   * time, NOT time steps.
   *
   * @param start time to start returning results at
   * @param stop time to stop returning results at
   * @param offset the amount of time to shift the signal being
   * applied to this excitation by. If the start time is non zero,
   * this can be used to shift the signal so that it ramps up
   * properly after the start time.
   */
  inline void set_time_param(float start, float stop, 
                             float offset)
  {
    time_start_ = start;
    time_stop_ = stop;
    time_offset_ = offset; 
  }

  /**
   * Returns the CSGBox the excitation is being applied to. Mostly
   * just for debugging.
   */ 
  inline shared_ptr<CSGBox> get_region()
  { return box_; }

  /**
   * Print a string representation to an ostream.
   */
  virtual ostream& to_string(ostream &os) const;

  friend ostream& operator<< (ostream& os, const Excitation &ex);
};

inline ostream& operator<< (ostream& os, const Excitation &ex)
{
  return ex.to_string(os);
}

#endif // EXCITE_H

