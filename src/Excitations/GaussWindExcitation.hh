/* 
   phred - Phred is a parallel finite difference time domain
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

#ifndef GAUSSWIND_EXCITATION_H
#define GAUSSWIND_EXCITATION_H

#include "WindowedExcitation.hh"

/**
 * An excitation that is windowed in 3 space by a Gaussian
 * window. This should make a pseudo plane wave...
 *
 * http://www.cg.tuwien.ac.at/studentwork/CESCG/CESCG99/TTheussl/node12.html
 */
class GaussWindExcitation : public WindowedExcitation
{
private:
  float std_dev_; /**< Standard deviation of Gaussian window */ 

protected:
  /**
   * Defines a Gaussian window function. 
   */
  virtual field_t window(region_t r, unsigned int i, unsigned int j, 
                         unsigned int k);

public:
  GaussWindExcitation(SourceFunction *sf); 

  ~GaussWindExcitation();
  
  void set_std_dev(float stddev);
  float get_std_dev();

};

#endif // GAUSSWIND_EXCITATION_H