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

#include <vector>
#include <iostream>
#include <string>

#include "Types.hh"
#include "Material.hh"

using namespace std; 

class MaterialLib {
 private:
  vector <Material> materials_;

 protected:

 public:
  // An input iterator class for traversing the materials list
  //class iterator {
  //private:
  //  vector<Material>::iterator iter_;
  //  
  //public:
  //  
  //};

  // Functions
  MaterialLib();
  ~MaterialLib();

  /**
   * Add a material to the material library. A copy of the material
   * object is stored, not a pointer or a reference.
   *
   * @param mat the material to store
   */
  void add_material(Material& mat);

  /**
   * Returns the number of materials in the library
   * @return an int
   */
  int num_materials();

  /**
   * Returns an interator that points to the start of the vector
   * containing the materials. TEMPORARY until a proper iterator for
   * the materials is set up.
   *
   * @return a vector<Material>::iterator
   */
  vector<Material>::iterator get_material_iter_begin();

  /**
   * Returns an interator that points to the end of the vector
   * containing the materials. TEMPORARY until a proper iterator for
   * the materials is set up.
   *
   * @return a vector<Material>::iterator
   */
  vector<Material>::iterator get_material_iter_end();
};

#endif // MATERIAL_LIB
