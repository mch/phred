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

#ifndef META_FDTD_H
#define META_FDTD_H

#include "FDTD.hh"

enum MetaType {
  METAFDTD_ONE, /**< One field component at a time */ 
  METAFDTD_THREE, /**< All E or all H field components at a time */ 
  METAFDTD_SGI_ORIGIN /**< Simple attempt at cache optimization on SGI
                         Origin machines, assuming 32kb L1 cache. */ 
  //METAFDTD_
};

/**
 * This is an FDTD class derived from the normal FDTD class. It uses
 * meta-programmable grid updates. Yay.
 */ 
class MetaFDTD : public FDTD
{
public:
  MetaFDTD();
  virtual ~MetaFDTD();

  /**
   * Run the simulation.
   */
  virtual void run();
    
  /**
   * Select the type of meta update to use
   */ 
  void set_metatype(MetaType mt)
  { mt_ = mt; }

  /**
   * Returns the type of meta update being used
   */ 
  MetaType get_metatype()
  { return mt_; }

private:
  /**
   * The method to use
   */
  MetaType mt_;

  region_t e_update_r_;
  region_t h_update_r_;

  /**
   * Apply the meta update to the E field
   */ 
  void update_e(GridUpdateData &gud);

  /**
   * Apply the meta update to the H field
   */ 
  void update_h(GridUpdateData &gud);

  /**
   * Computes the E and H regions where all three field components can
   * be computed at one time.
   */ 
  void compute_update_regions();
};

#endif // META_FDTD_H
