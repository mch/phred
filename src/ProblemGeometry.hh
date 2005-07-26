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

#ifndef PROBLEM_GEOMETRY_H
#define PROBLEM_GEOMETRY_H

#include "LifeCycle.hh"
#include "CSG/CSGObject.hh"
#include "CSG/CSGBox.hh"

#include <string>
#include <vector>
using namespace std;

#include <sys/types.h>
#include <boost/shared_ptr.hpp>
using namespace boost;

/**
 * "Metadata" attached to CSGObjects
 */ 
class GeomObject {
public:
  GeomObject(string material, shared_ptr<CSGObject> obj)
    : obj_(obj), material_(material)
  {}

  shared_ptr<CSGObject> obj_;
  string material_;
  unsigned int material_id_;
};

/**
 * This class is a container that holds information about the
 * geometric properties of this problem. It has box object which
 * represents the grid, and any number of other CSG objects which
 * make up the features within the grid. 
 */
class ProblemGeometry : public LifeCycle {
public:
  ProblemGeometry();
  ~ProblemGeometry();

  /**
   * Set the grid material
   */ 
  void set_grid_material(const char *material);

  /**
   * Get the default material index inside the grid. 
   */ 
  mat_idx_t get_grid_material_id() const;

  /**
   * Returns the material id for a given point in the grid. 
   */ 
  unsigned int get_material_id(float x, float y, float z) const;

  /**
   * Add a CSG object
   */ 
  void add_object(string material, shared_ptr<CSGObject> obj);

  /**
   * Set the grid size
   */
  void set_grid_size(float x_size, float y_size, float z_size);

  /**
   * Returns the grid size
   */ 
  point get_grid_size() const;

  /**
   * Set the position of the centre of the grid. 
   */ 
  void set_grid_centre(float x, float y, float z);

  /**
   * Returns the position of the centre of the grid
   */ 
  point get_grid_centre() const;

  /**
   * Set the units multiplier. The user can choose to rescale the
   * objects, or let them remain the same size. I.e. If the old units
   * is meters, and it is being changed to mm, scaling will make a 1m
   * diameter sphere into a 1mm diameter sphere. 
   */
  //void set_units(float multiplier, bool scale);

  /**
   * Initialize the problem geometry; get the material id numbers for
   * the geometries. 
   */
  void init(const Grid &grid);

  /**
   * Deinit
   */ 
  void deinit();

protected:
  vector <GeomObject> objects_;

  // Grid box
  GeomObject grid_box_;

  // Units multiplier
  float unit_;
};

#endif // PROBLEM_GEOMETRY_H
