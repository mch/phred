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

#include "WaveguideExcitation.hh"

WaveGuideExcitation::WaveGuideExcitation(SouceFunction *sf)
  : WindowedExcitation(sf)
{}

WaveGuideExcitation::~WaveGuideExcitation()
{}

field_t WaveGuideExcitation::window(region_t r, 
                                    unsigned int x, 
                                    unsigned int y, 
                                    unsigned int z)
{
  
}

