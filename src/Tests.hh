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

/**
 * \file Tests.hh
 * Compiled in tests for various situations. Will probably be
 * eliminated once Jan's script loader is complete. 
 */

// void point_test(int rank, int size);
// void pml_test(int rank, int size);
// void takakura_test(int rank, int size);
// void laser_test(int rank, int size);
// void coupler_test(int rank, int size);



void hole();
void mn_benchmark();
void var_benchmark(unsigned int x_cells, unsigned int y_cells, 
                   unsigned int z_cells);
void grooves_top();
void grooves_bottom();
void grooves_both();

void square_hole(int ysize);
void square_hole_Ag(int ysize);
