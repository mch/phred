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

#ifndef INTERVAL_H
#define INTERVAL_H

#include <ostream>

/**
 * This class defines an closed interval of equally spaced
 * points. Primary use is in Results, for calculating a range of
 * frequencies or angles.
 *
 * \bug needs a const iterator... Use boost's iterator_facade:
 * http://boost.org/libs/iterator/doc/iterator_facade.html
 */ 
template<class T>
class Interval 
{
public:
  /**
   * Default constructor
   */ 
  Interval()
    : start_(0), end_(0), space_(0), items_(0), num_pts_(0)
  {}

  /**
   * Destructor
   */ 
  ~Interval()
  {
    if (items_)
      delete[] items_;
  }

  /**
   * Copy constructor
   */ 
  Interval(const Interval<T> &rhs)
    : items_(0)
  {
    start_ = rhs.start_;
    end_ = rhs.end_;
    num_pts_ = rhs.num_pts_;

    calc_items();
  }

  /**
   * Assignment operator
   */ 
  Interval<T> &operator=(const Interval<T> &rhs)
  {
    start_ = rhs.start_;
    end_ = rhs.end_;
    num_pts_ = rhs.num_pts_;

    calc_items();

    return *this;
  }

  /**
   * Constructor which initializes the object.
   *
   * @param start starting point
   * @param stop end point
   * @param num number of points in the interval. Must be > 0.
   */ 
  Interval(const T &start, const T &end, unsigned int &num)
    : start_(start), end_(end), items_(0), num_pts_(num)
  {
    calc_items();
  }

  /**
   * Set the interval parameters
   *
   * @param start starting point
   * @param end end point
   * @param num number of points in the interval. Must be > 0.
   */ 
  void set_params(const T &start, const T &end, unsigned int &num)
  {
    start_ = start;
    end_ = end;
    num_pts_ = num;

    calc_items();
  }

  /**
   * Get the starting point of the interval.
   */ 
  inline const T &get_start() const
  { return start_; }

  /**
   * Get the end point of the interval.
   */ 
  inline const T &get_end() const
  { return end_; }
  
  /**
   * Get the number of items in the interval.
   */ 
  inline const unsigned int &length() const
  { return num_pts_; }

  /**
   * Get the spacing between two points. 
   */ 
  inline const T &get_spacing() const
  { return space_; }

  /**
   * Get an item in the interval. 
   *
   * @param item the number of the item to get. Must be >= zero
   * and < the value returned by length(). If it is not in this
   * interval, 0 will be returned.
   */ 
  inline T get(unsigned int item)
  { 
    if (item >= num_pts_ || !items_)
      return 0;

    return items_[item];
  }

  /**
   * Returns a pointer to the array of items making up this
   * Interval. BE VERY CAREFUL WITH THIS POINTER! THIS BREAKS
   * ENCAPSULATION!
   */ 
  inline const T *get_ptr()
  { return items_; }

protected:
  T start_;              /**< starting point */ 
  T end_;                /**< end point */ 
  T space_;              /**< spacing between two items */ 
  T *items_;             /**< items in the interval */ 
  unsigned int num_pts_; /**< Number of points in the interval */ 

  /**
   * Calculate the items in the interval. This function checks the
   * parameters for validity before building the list of items. If it
   * is flawed, then the length of the interval is set to zero.
   */ 
  void calc_items()
  {
    if (end_ < start_ || num_pts_ < 1)
    {
      num_pts_ = 0;
      return;
    }
   
    if (num_pts_ > 0)
      space_ = (end_ - start_) / (num_pts_ - 1);
    else 
      space_ = 0;

    if (items_)
      delete[] items_;

    items_ = new T[num_pts_];

    for (unsigned int i = 0; i < num_pts_; i++)
      items_[i] = start_ + i * space_;
  }

  template<class S>
  friend std::ostream &operator<<(std::ostream &os, 
                           const Interval<S> &interval);
};

template<class T>
std::ostream &operator<<(std::ostream &os, const Interval<T> &interval)
{
  os << "[";

  if (num_pts_ > 0 && items_)
  {
    for (unsigned int i = 0; i < num_pts_ - 1; i++)
    {
      os << items_[i] << " ";
    }

    os << items_[num_pts_ - 1];
  }

  os << "]";

  return os;
}

#endif // INTERVAL_H
