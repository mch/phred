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

/** \class MaterialLib
 * \brief A material library.
 * This class is a material library for keeping track of materials. 
 * It contains a vector of material objects, and is able to generate
 * arrays of material coefficients for use in the FDTD algorithm. 
 */

#ifndef MATERIAL_LIB
#define MATERIAL_LIB

#include "config.h"

#include <vector>
#include <map>
#include <iostream>
#include <string>

#include "Types.hh"
#include "Material.hh"

using namespace std; 

/**
 * This class represents a collection of materials. Materials
 * representing perfect electric conductor (PEC) and free space
 * (freespace) are present in a newly created library by default. 
 *
 * It can load and save material information from and to a text
 * file. The format for the text file is as follows:
 *
 * The first line is a header line descibing the file. This is
 * followed by an empty line, and the Material definition
 * blocks, which are seperated by an empty line. Each block consists
 * of a string naming the material, then
 * a space, then four floating point numbers identifying the
 * permittivity, permeability, conductivity, and magnetic
 * conductivity of the medium. 
 *
 * Materials may have additional properties. These are specified by
 * additional lines following the material block. The first item on
 * the line is a string identifying the property, and the second is
 * the floating point value of the property. 
 *
 * An example file:
 *
 * Phred 0.1.0 Material Library
 * 
 * glass 1.8 1 0 0
 *
 * gold 1 1 0 0
 * plasma_freq 2.4e-12
 * collision_freq 2.4e-12
 *
 * silver 1 1 0 0
 * plasma_freq 2e15
 * collision_freq 5.7e13
 *
 *
 *
 * The above collision frequency and plasma frequencies must be set in
 * rad/sec.
 */

class MaterialLib {
  friend class Grid;
  friend class UPmlCommon;

 private:
  map <string, Material> materials_;

 protected:

 public:
  MaterialLib();
  ~MaterialLib();

  /**
   * Add a material to the material library. A copy of the material
   * object is stored, not a pointer or a reference.
   *
   * @param mat the material to store
   */
  void add_material(const char *name, Material& mat);

  /**
   * Return a reference to a material object
   */
  const Material &get_material(const char *name) const;

  /**
   * Returns the number of materials in the library
   * @return an int
   */
  int num_materials() const;

  // temp
  inline map<string, Material>::const_iterator get_material_iter_begin() const
  { return materials_.begin(); }

  inline map<string, Material>::const_iterator get_material_iter_end() const
  { return materials_.end(); }


  /**
   * Load material information from a text file. 
   */ 
  void load_material_file(const char *filename);

  /**
   * Save the material information to a text file. 
   */ 
  void save_material_file(const char *filename);

};

#endif // MATERIAL_LIB
