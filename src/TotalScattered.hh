/* 
   phred - Phred is a parallel finite difference time domain
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

#ifndef TOTAL_SCATTERED_H
#define TOTAL_SCATTERED_H

/**
 * This class implements the total-scattered field forumlation. It
 * defines a region of the FDTD grid as the total field area, where
 * the field components inside the grid are defined to be the total
 * field. Outside this region, the incident field is subtracted and
 * the field components are considered to consist of only the
 * scattered fields. 
 */
class TotalScattered 
{
private:
protected:
public:
  TotalScattered(region_t region);
  ~TotalScattered();

  
};

#endif
