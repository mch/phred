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

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "LifeCycle.hh"
#include "Types.hh"

/**
 * This is an abstract base class which defines an interface which is
 * to be used to build objects which define... geometry
 * objects. Boxes and spheres and stuff. You know. Each object must
 * know how to iterate over the cells it constrains, so that you can
 * do things to it. 
 *
 * Geometry objects can store auxilary data that may be required for
 * dispersion relationships. 
 */
class Geometry : public LifeCycle
{
private:
protected:
  unsigned int material_id_;

  field_t **aux_; /**< Auxilariy data required for dispersion
                     relationships that happen inside the geometry. */ 

  region_t bounding_box_; /**< A box that entierly contains the object
                             defined by any given geometry. */  

  region_t local_bb_; /**< Bounding box in local coordinates */ 

public:
  Geometry() {}
  virtual ~Geometry() {}

  /**
   * Iterator class for traversing the cells within the region defined
   * by the geometry. 
   */
  class iterator {
  private:
  protected:
  public:
    iterator();
    ~iterator();

    /**
     * Postfix iterator increment
     */
    const iterator &operator++();

    /**
     * Prefix iterator increment
     */
    const iterator &operator++(int a); // , const iterator &rhs);
    
    /**
     * Dereference the iterator; return a index into the grid where
     * the current element can be found. 
     */
    unsigned int operator*(const iterator &rhs);
  };

  /**
   * Set the auxiliary field data
   */
  inline void set_aux_data(field_t **aux)
  {
    aux_ = aux;
  }

  /** 
   * Return the aux field data
   */
  inline field_t **get_aux_data()
  {
    return aux_;
  }

  /**
   * Set the material id to use
   */
  inline void set_material_id(unsigned int mid)
  {
    material_id_ = mid;
  }

  /**
   * Returns the material ID for this geometry
   */
  inline unsigned int get_material_id()
  {
    return material_id_;
  }

  /**
   * Init the geometry against the grid; check that the geometry is
   * inside the grid
   */ 
  virtual void init(const Grid &grid);

  /**
   * Deinit the geometry against the grid
   */ 
  virtual void deinit();

  /**
   * Update the material indicies of the grid 
   */
  virtual void set_material(Grid &grid) = 0;

  /**
   * Returns true if the given point is inside this geometry, with
   * respect to local grid coordinates. 
   */ 
  virtual bool local_point_inside(unsigned int x,
                                  unsigned int y, 
                                  unsigned int z) = 0;

  /**
   * Returns true if the given point is inside this geometry, with
   * respect to global grid coordinates. 
   */ 
  //virtual bool global_point_inside(unsigned int x,
  //                                 unsigned int y, 
  //                                 unsigned int y) = 0;

  /**
   * Returns the region which defines the bounding box for this object.
   */
  inline region_t get_bounding_box()
  {
    return bounding_box_;
  }

  /**
   * Returns the region which defines the bounding box for this object
   * in local coordinates.
   */
  inline region_t get_local_bounding_box()
  {
    return local_bb_;
  }

};

#endif // GEOMETRY_H
