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

#ifndef CONTIGUOUS_H
#define CONTIGUOUS_H

/**
 * This class contains functions for dealing with a contiguous chunk
 * of memory which represents an N dimensional array of data.
 *
 * \bug This needs to be smarter when dealing with contiguous
 * memory... it should detect contiguous blocks and avoid creating a
 * huge number of types.
 */ 
//template<class T, int N = 3>
class Contiguous {
public:

  /** 
   * Calculate the index of a point at coord in a contiguous array
   * representing an N dimensional matrix of size size. Assumes data is
   * contiguous along the last dimension.
   *
   * @param N the number of dimensions
   * @param size the size of each dimension
   * @param coord the point to calculate the index for. 
   */ 
  static unsigned int idx(const unsigned int N, const unsigned int *size, 
                          const unsigned int *coord)
  {
    unsigned int idx = coord[0];

    for (int i = 0; i < N - 1; i++)
      idx = coord[i + 1] + idx * size[i + 1];
    
    return idx;
  }

  /**
   * Compute displacements for an indexed data type that pulls out a
   * chunk of a larger N dimensional matrix. This assumes that the data
   * varies fastest (is contiguous along) the last dimension specified.
   *
   * @param N The number of dimensions in the array
   * @param dsize The size of each dimension
   * @param coord The starting coordinate
   * @param lens The length each dimension must span
   * @return num_displacements A pointer to an integer where this
   * function can store the number of displacements that have been
   * computed.
   * @return displacements The pointer to a pointer which this
   * function will set to the address of a new array of ints
   * containing the calculated displacements. The client MUST free
   * this memory. 
   */ 
  static void compute_displacements(const unsigned int N, 
                                    const unsigned int *dsize, 
                                    const unsigned int *coord, 
                                    const unsigned int *lens,
                                    unsigned int *num_displacements,
                                    int **displacements)
  {
    unsigned int num_disps = 1;
    unsigned int disp_idx = 1;
    unsigned int offset = dsize[N - 1];
    unsigned int negoffset = 1;
    
    for (int i = 0; i < N - 1; i++)
    {
      num_disps *= lens[i];
    }
    
    int *disps = new int [num_disps];
    
    disps[0] = idx(N, dsize, coord);
    
    disp_idx = displ_helper(N, dsize, lens, 1, disp_idx, disps, &offset, 
                            &negoffset, true);

    *displacements = disps;
    *num_displacements = num_disps;
  }

private:
  /**
   * Displacement computation recursion helper
   */ 
  static int displ_helper(const unsigned int N, 
                          const unsigned int *dsize, 
                          const unsigned int *lens,
                          const unsigned int curr_dim, 
                          unsigned int disp_idx, 
                          int *disps, 
                          unsigned int *offset,
                          unsigned int *negoffset, 
                          const bool first)
  {
    // Base case
    if (curr_dim == N)
    {
      return disp_idx;
    }
    
    int j = 0;
    
    if (first) {
      j = 1;

      disp_idx = displ_helper(N, dsize, lens, curr_dim + 1, disp_idx, disps, 
                              offset, negoffset, true);
    }

    for (; j < lens[curr_dim - 1]; j++)
    {
      if (curr_dim == N - 1)
      {
        disps[disp_idx] = disps[disp_idx - *negoffset] + *offset;
        disp_idx++;
      } else {
        disp_idx = displ_helper(N, dsize, lens, curr_dim + 1, disp_idx, disps, 
                                offset, negoffset, false);
      }
    }
    
    if (first)
    {
      *negoffset = lens[curr_dim - 1];
      *offset *= dsize[curr_dim - 1];
    }

    return disp_idx;
  }                          


};

#endif // CONTIGUOUS_H
