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

#ifndef BARTLETT_EXCITATION_H
#define BARTLETT_EXCITATION_H

#include "WindowedExcitation.hh"

/**
 * An excitation that is windowed in 3 space by Bartlett windows (tent
 * functions). This should make a pseudo plane wave...
 * 
 * Check the Bartlett Window definition at:
 * http://www.cg.tuwien.ac.at/studentwork/CESCG/CESCG99/TTheussl/node6.html
 */
class BartlettExcitation : public WindowedExcitation
{
private:
protected:
  /**
   * Defines a bartlett window function. 
   */
  virtual field_t window(float x, float y, float z);

public:
  BartlettExcitation(SourceFunction *sf);
  ~BartlettExcitation();
  
};

#endif // BARTLETT_EXCITATION_H
