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

#ifndef PHRED_REGION_H
#define PHRED_REGION_H

#include "Exceptions.hh"
#include "GridInfo.hh"

/**
 * A region represents a block of Yee cells in a FDTD computational
 * domain. Since a single computational domain may be de-composed into
 * multiple domains which may then be assigned to different nodes or
 * processors, it is necessary to track both global and local
 * information about a region. The local part of the region refers to
 * the block of Yee cells that are actually allocated and located on
 * the current processor.
 * 
 * A region is a seperate concept from a geometry (CSG) object such as
 * a box, in that geometry objects are given in terms of real numbers
 * which represent a physical distance, while regions are unsigned
 * integers which refer to Yee cells.
 *
 * Regions specifed to be *inclusive*; xmin is the index of the Yee
 * cell closest to the origin that is included in the region, and xmax
 * is the index of the Yee cell farthest from the origin. This is
 * noted because it is changed from older versions of Phred and from
 * the grandfather C code. 
 *
 * This region NEVER includes overlapping Yee cells, which arise when
 * the grid is split across nodes. 
 *
 * The global part of the region is immutable. Min and max values are
 * ALWAYS set in global coordinates. 
 */
class Region {
public:
  /**
   * Class contructor. The region starts at 0,0,0 and contains 0 cells
   * along each dimension. The region does not have a parent; global
   * and local values are the same.
   */
  Region();

  /**
   * Class contructor. The region starts at 0,0,0 and contains 0 cells
   * along each dimension. The given global region is interperted as
   * the region containing this one; it is a superset of this
   * region. A copy of the global is stored.
   */
  Region(const Region &global);

  /**
   * Class contructor. This one initializes the region to a non-null space. 
   * The region does not have a parent; global and local values are
   * the same.
   */
  Region(unsigned int xmin, unsigned int xmax,
         unsigned int ymin, unsigned int ymax,
         unsigned int zmin, unsigned int zmax);

  /**
   * Class contructor. This one initializes the region to a non-null space. 
   * The region does not have a parent; global and local values are
   * the same. The given global region is interperted as
   * the region containing this one; it is a superset of this
   * region. A copy of the global is stored.
   */
  Region(unsigned int xmin, unsigned int xmax,
         unsigned int ymin, unsigned int ymax,
         unsigned int zmin, unsigned int zmax,
         const Region &global);

  virtual ~Region();

  /**
   * Indicates if this region has any cells in the grid that is local
   * to this processor. This should be called before calling any of
   * the get_{x,y,z}{min,max}() functions. 
   */ 
  virtual bool is_local() const;

  /**
   * Set the closed interval along the x axis containing this region. 
   *
   * @param xmin
   * @param xmax
   * @exception ResultException if xmin > xmax. 
   */ 
  void set_x(unsigned int xmin, unsigned int xmax);

  /**
   * Set the closed interval along the y axis containing this region. 
   *
   * @param ymin
   * @param ymax
   * @exception ResultException if ymin > ymax. 
   */ 
  void set_y(unsigned int ymin, unsigned int ymax);

  /**
   * Set the closed interval along the z axis containing this region. 
   *
   * @param zmin
   * @param zmax
   * @exception ResultException if zmin > zmax. 
   */ 
  void set_z(unsigned int zmin, unsigned int zmax);

  /**
   * Returns the global xmin of this region. 
   */
  inline unsigned int get_global_xmin() const
  { return global_xmin_; }

  /**
   * Returns the global ymin of this region. 
   */
  inline unsigned int get_global_ymin() const
  { return global_ymin_; }

  /**
   * Returns the global zmin of this region. 
   */
  inline unsigned int get_global_zmin() const
  { return global_zmin_; }

  /**
   * Returns the global xmax of this region. 
   */
  inline unsigned int get_global_xmax() const
  { return global_xmax_; }

  /**
   * Returns the global ymax of this region. 
   */
  inline unsigned int get_global_ymax() const
  { return global_ymax_; }

  /**
   * Returns the global zmax of this region. 
   */
  inline unsigned int get_global_zmax() const
  { return global_zmax_; }
  



  /**
   * Returns the local xmin of this region. This is *VIRTUAL* and NOT
   * inlined, so for high performance code like loops, store this
   * result in a variable.
   */
  virtual unsigned int get_xmin() const;

  /**
   * Returns the local ymin of this region. This is *VIRTUAL* and NOT
   * inlined, so for high performance code like loops, store this
   * result in a variable.
   */
  virtual unsigned int get_ymin() const;

  /**
   * Returns the local zmin of this region. This is *VIRTUAL* and NOT
   * inlined, so for high performance code like loops, store this
   * result in a variable.
   */
  virtual unsigned int get_zmin() const;

  /**
   * Returns the local xmax of this region. This is *VIRTUAL* and NOT
   * inlined, so for high performance code like loops, store this
   * result in a variable.
   */
  virtual unsigned int get_xmax() const;

  /**
   * Returns the local ymax of this region. This is *VIRTUAL* and NOT
   * inlined, so for high performance code like loops, store this
   * result in a variable.
   */
  virtual unsigned int get_ymax() const;

  /**
   * Returns the local zmax of this region. This is *VIRTUAL* and NOT
   * inlined, so for high performance code like loops, store this
   * result in a variable.
   */
  virtual unsigned int get_zmax() const;
  
protected:
  unsigned int global_xmin_;
  unsigned int global_xmax_;

  unsigned int global_ymin_;
  unsigned int global_ymax_;

  unsigned int global_zmin_;
  unsigned int global_zmax_;

  Region *global_;
};

/**
 * This region class INCLUDES the any subdomain overlap that arises
 * when a a Grid is split across nodes. 
 */ 
class OverlapRegion : public Region
{
public:
  /**
   * Class contructor. 
   */
  OverlapRegion(const Region &global);

  /**
   * Class contructor. This one initialized the region to a non-null space.  
   */
  OverlapRegion(unsigned int xmin, unsigned int xmax,
                unsigned int ymin, unsigned int ymax,
                unsigned int zmin, unsigned int zmax, 
                const Region &global);

  virtual ~OverlapRegion();

  /**
   * Indicates if this region has any cells in the grid that is local
   * to this processor. This should be called before calling any of
   * the get_{x,y,z}{min,max}() functions. 
   */ 
  virtual bool is_local() const;

   /**
   * Returns the local xmin of this region. This is *VIRTUAL* and NOT
   * inlined, so for high performance code like loops, store this
   * result in a variable.
   */
  virtual unsigned int get_xmin() const;

  /**
   * Returns the local ymin of this region. This is *VIRTUAL* and NOT
   * inlined, so for high performance code like loops, store this
   * result in a variable.
   */
  virtual unsigned int get_ymin() const;

  /**
   * Returns the local zmin of this region. This is *VIRTUAL* and NOT
   * inlined, so for high performance code like loops, store this
   * result in a variable.
   */
  virtual unsigned int get_zmin() const;

  /**
   * Returns the local xmax of this region. This is *VIRTUAL* and NOT
   * inlined, so for high performance code like loops, store this
   * result in a variable.
   */
  virtual unsigned int get_xmax() const;

  /**
   * Returns the local ymax of this region. This is *VIRTUAL* and NOT
   * inlined, so for high performance code like loops, store this
   * result in a variable.
   */
  virtual unsigned int get_ymax() const;

  /**
   * Returns the local zmax of this region. This is *VIRTUAL* and NOT
   * inlined, so for high performance code like loops, store this
   * result in a variable.
   */
  virtual unsigned int get_zmax() const;

protected:

private:
};

#endif // PHRED_REGION_H
