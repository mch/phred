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

#ifndef JAN_FDTD_H
#define JAN_FDTD_H

#include "FDTD.hh"

/**
 * This class is used by the parse for Jan's grammer to store the
 * objects it generates. The only real difference between this class
 * and FDTD is that this class takes responsibility for cleaning up
 * the memory of the objects it knows about, and a bunch of functions
 * are provided that do input error checking. 
 *
 * A modified version of Jan's grammer is used which stores everything
 * about the problem in one file and takes the relationship between
 * Results and DataWriters implemented in this program into account. 
 *
 * \section jangrammer Parsing Jan's Grammer
 * The input file that is accepted by phred is a modified version of
 * Jan's grammer. It is parsed by flex and objects are created by a
 * yacc
 */ 
class JanFDTD : public FDTD
{
private:
protected:
public:
  JanFDTD();
  ~JanFDTD();

  virtual void parse_file(string filename);

  void set_program_mode(char m);
  void set_structure_mode(char m);
  void set_timestep_mode(char m);
  void set_runtime(unsigned int rt);
  void set_time_modulo(unsigned int rt);

  void set_dimx(unsigned int dx);
  void set_dimy(unsigned int dy);
  void set_dimz(unsigned int dz);

  void set_deltax(double dx);
  void set_deltay(double dy);
  void set_deltaz(double dz);
  void set_deltat(double dt);

  void set_num_materials(unsigned int nm);
  void add_material(unsigned int id, float eps, float sigma, 
                    float mu, float sigma_star);

  void set_pml_boundary(Face face, unsigned int thickness, 
                        char p, unsigned int nrml_refl);

  void set_boundary(Face face, BoundaryCondition bc);
};

#endif // JAN_FDTD_H
