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

#ifndef WINDOWED_EX_RESULT_H
#define WINDOWED_EX_RESULT_H

#include "Result.hh"
#include "../Excitations/WindowedExcitation.hh"

/**
 * This class outputs the spatial window which is applied to a signal
 * by a subclass of WindowedExcitation.
 *
 * This is pretty much just for testing. 
 */
class WindowedExResult : public Result
{
public:
  WindowedExResult(shared_ptr<WindowedExcitation> wex);

  void calculate_pre_result(const Grid &grid);

  void init(const Grid &grid);
  void deinit();
  
protected:
  float *wnd_;

  Variable var_;
  shared_ptr<WindowedExcitation> wex_;
  shared_ptr<Block> gregion_;
};

#endif // WINDOWED_EX_RESULT_H
