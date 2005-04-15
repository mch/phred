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

#ifndef PROP_RESULT_H
#define PROP_RESULT_H

#include "Result.hh"

/**
 * This class returns things like grid and time discritization
 * information and other small bits of information that may be useful
 * in post-processing.
 */ 
class PropertiesResult : public Result
{
public:
  /**
   * Properties are made available here. 
   */ 
  void calculate_pre_result(const Grid &grid);

  /**
   * Initalize the result
   */
  void init(const Grid &grid);

  /**
   * Deinit the result
   */ 
  void deinit();

  /**
   * Print a string representation to an ostream.
   */
  ostream& to_string(ostream &os) const;

private:
  Variable dx_var_; /**< Grid x delta */  
  Variable dy_var_; /**< Grid y delta */  
  Variable dz_var_; /**< Grid z delta */  
  Variable dt_var_; /**< Time delta */  

  Variable ts_var_; /**< Number of time steps */  

  delta_t dx_, dy_, dz_, dt_;

  unsigned int dimx_, dimy_, dimz_, ts_;
};

#endif // PROP_RESULT_H
