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

#ifndef GRID_PLANE_H
#define GRID_PLANE_H

#include "Grid.hh"

/**
 * This is an abstract base class which serves as a template adapter
 * class for accessing normal and tangential field components on grid
 * planes. All functions are inline, so they should be optimized
 * out. The use of these adapter classes serves to eliminate the need
 * to return pointers to protected members of the Grid class, which
 * is bad OO karma. This keeps our OO karma points high, and ensures
 * that the assert statements used in the set/get methods of the grid
 * class are honored in DEBUG mode. That really the most important
 * part; I'd rather see an assert message than an intermittant
 * segfault and garbage data. 
 *
 * Adapter classes can be used polymorphically inside boundary
 * condition loops. If polymorphism is too slow, they can be passed
 * as template paramters to a templated loop. 
 *
 * Components are related thusly: the cross product of the first and
 * second tangential components points in the direction of normal
 * component. For counting tangential components, it's X, Y, then
 * Z. If Z is first, the second is X.
 */
class GridPlane
{
private:
protected:
  Grid &grid_; /**< a reference to the grid... performance implications? */

public:
  /**
   * Constructor; the grid to Adapt must be passed in.
   */
  GridPlane(Grid &grid)
    : grid_(grid)
  {}

  /**
   * Destructor
   */
  virtual ~GridPlane()
  {}

  /**
   * Returns a pointer the normal E component. This is intended to be
   * used in Result gathering and update loops to make things faster.
   * 
   * @param x
   * @param y
   * @param z
   */ 
  virtual const field_t *get_e_n_ptr(unsigned int x, unsigned int y, 
                                     unsigned int z) const = 0;

  /**
   * Returns a pointer the normal H component. This is intended to be
   * used in Result gathering and update loops to make things faster.
   * 
   * @param x
   * @param y
   * @param z
   */ 
  virtual const field_t *get_h_n_ptr(unsigned int x, unsigned int y, 
                                     unsigned int z) const = 0;

  /**
   * Returns a pointer the first tangential E component. This is
   * intended to be used in Result gathering and update loops to make
   * things faster.
   * 
   * @param x
   * @param y
   * @param z
   */ 
  virtual const field_t *get_e_t1_ptr(unsigned int x, unsigned int y, 
                                      unsigned int z) const = 0;

  /**
   * Returns a pointer the first tangential H component. This is
   * intended to be used in Result gathering and update loops to make
   * things faster.
   * 
   * @param x
   * @param y
   * @param z
   */ 
  virtual const field_t *get_h_t1_ptr(unsigned int x, unsigned int y, 
                                      unsigned int z) const = 0;

  /**
   * Returns a pointer the second tangential E component. This is
   * intended to be used in Result gathering and update loops to make
   * things faster.
   * 
   * @param x
   * @param y
   * @param z
   */ 
  virtual const field_t *get_e_t2_ptr(unsigned int x, unsigned int y, 
                                      unsigned int z) const = 0;

  /**
   * Returns a pointer the second tangential H component. This is
   * intended to be used in Result gathering and update loops to make
   * things faster.
   * 
   * @param x
   * @param y
   * @param z
   */ 
  virtual const field_t *get_h_t2_ptr(unsigned int x, unsigned int y, 
                                      unsigned int z) const = 0;

  /** 
   * Set the normal E field component
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  virtual void set_e_n(unsigned int x, unsigned int y, 
                       unsigned int z, field_t value) = 0;

  /** 
   * Set the normal H field component
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  virtual void set_h_n(unsigned int x, unsigned int y, 
                       unsigned int z, field_t value) = 0;

  /**
   * Set the first tangential E field component, E_1. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  virtual void set_e_t1(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value) = 0;

  /**
   * Set the second tangential E field component, E_2. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  virtual void set_e_t2(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value) = 0;

  /**
   * Set the first tangential H field component, H_1. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  virtual void set_h_t1(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value) = 0;

  /**
   * Set the second tangential H field component, H_2. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  virtual void set_h_t2(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value) = 0;


  /** 
   * Get the normal E field component
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  virtual field_t get_e_n(unsigned int x, unsigned int y, 
                          unsigned int z) const = 0;

  /** 
   * Get the normal H field component
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  virtual field_t get_h_n(unsigned int x, unsigned int y, 
                          unsigned int z) const = 0;

  /**
   * Get the first tangential E field component, E_1. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  virtual field_t get_e_t1(unsigned int x, unsigned int y, 
                           unsigned int z) const = 0;

  /**
   * Get the second tangential E field component, E_2. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  virtual field_t get_e_t2(unsigned int x, unsigned int y, 
                           unsigned int z) const = 0;

  /**
   * Get the first tangential H field component, H_1. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  virtual field_t get_h_t1(unsigned int x, unsigned int y, 
                           unsigned int z) const = 0;

  /**
   * Get the second tangential H field component, H_2. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  virtual field_t get_h_t2(unsigned int x, unsigned int y, 
                           unsigned int z) const = 0;


  /**
   * Returns an average of 2 Et1 field components so that the result
   * is the value of the returned field is what would be present at the
   * centre of the cell face. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return averaged field component value
   */
  virtual field_t get_avg_e_t1(unsigned int x, unsigned int y, 
                               unsigned int z) const = 0;
   
  /**
   * Returns an average of 2 Et2 field components so that the result
   * is the value of the returned field is what would be present at the
   * centre of the cell face. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return averaged field component value
   */
  virtual field_t get_avg_e_t2(unsigned int x, unsigned int y, 
                               unsigned int z) const = 0;

  /**
   * Returns an average of 4 Ht1 field components so that the result
   * is the value of the returned field is what would be present at the
   * centre of the cell face. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return averaged field component value
   */
  virtual field_t get_avg_h_t1(unsigned int x, unsigned int y, 
                               unsigned int z) const = 0;
   
  /**
   * Returns an average of 4 Ht2 field components so that the result
   * is the value of the returned field is what would be present at the
   * centre of the cell face. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return averaged field component value
   */
  virtual field_t get_avg_h_t2(unsigned int x, unsigned int y, 
                               unsigned int z) const = 0;
   


};

/**
 * Grid adapter class for field components on the X = constant or YZ
 * plane. 
 *
 * Et1 = EY
 * Et2 = EZ
 * Ht1 = HY
 * Ht2 = HZ
 */
class YZPlane : public GridPlane
{
private:
protected:
public:
  YZPlane(Grid &grid)
    : GridPlane(grid)
  {}

  ~YZPlane()
  {}
  
  inline const field_t *get_e_n_ptr(unsigned int x, unsigned int y, 
                                    unsigned int z) const
  { return grid_.get_pointer(grid_point(x, y, z), FC_EX); }

  inline const field_t *get_h_n_ptr(unsigned int x, unsigned int y, 
                                    unsigned int z) const
  { return grid_.get_pointer(grid_point(x, y, z), FC_HX); }

  inline const field_t *get_e_t1_ptr(unsigned int x, unsigned int y, 
                                     unsigned int z) const
  { return grid_.get_pointer(grid_point(x, y, z), FC_EY); }

  inline const field_t *get_h_t1_ptr(unsigned int x, unsigned int y, 
                                     unsigned int z) const
  { return grid_.get_pointer(grid_point(x, y, z), FC_HY); }

  inline const field_t *get_e_t2_ptr(unsigned int x, unsigned int y, 
                                     unsigned int z) const
  { return grid_.get_pointer(grid_point(x, y, z), FC_EZ); }

  inline const field_t *get_h_t2_ptr(unsigned int x, unsigned int y, 
                                     unsigned int z) const
  { return grid_.get_pointer(grid_point(x, y, z), FC_HZ); }

  inline void set_e_n(unsigned int x, unsigned int y, 
                       unsigned int z, field_t value)
  { grid_.set_ex(x, y, z, value); }

  inline void set_h_n(unsigned int x, unsigned int y, 
                       unsigned int z, field_t value)
  { grid_.set_hx(x, y, z, value); }

  inline void set_e_t1(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  { grid_.set_ey(x, y, z, value); }

  inline void set_e_t2(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  { grid_.set_ez(x, y, z, value); }

  inline void set_h_t1(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  { grid_.set_hy(x, y, z, value); }

  inline void set_h_t2(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  { grid_.set_hz(x, y, z, value); }

  inline field_t get_e_n(unsigned int x, unsigned int y, 
                          unsigned int z) const
  { return grid_.get_ex(x, y, z); }

  inline field_t get_h_n(unsigned int x, unsigned int y, 
                          unsigned int z) const
  { return grid_.get_hx(x, y, z); }

  inline field_t get_e_t1(unsigned int x, unsigned int y, 
                           unsigned int z) const
  { return grid_.get_ey(x, y, z); }

  inline field_t get_e_t2(unsigned int x, unsigned int y, 
                           unsigned int z) const
  { return grid_.get_ez(x, y, z); }

  inline field_t get_h_t1(unsigned int x, unsigned int y, 
                           unsigned int z) const
  { return grid_.get_hy(x, y, z); }

  inline field_t get_h_t2(unsigned int x, unsigned int y, 
                           unsigned int z) const
  { return grid_.get_hz(x, y, z); }

  inline field_t get_avg_e_t1(unsigned int x, unsigned int y, 
                               unsigned int z) const
  { return 0.5 * (grid_.get_ey(x, y, z) + grid_.get_ey(x, y, z + 1)); }

  inline field_t get_avg_e_t2(unsigned int x, unsigned int y, 
                               unsigned int z) const
  { return 0.5 * (grid_.get_ez(x, y, z) + grid_.get_ez(x, y + 1, z)); }

  inline field_t get_avg_h_t1(unsigned int x, unsigned int y, 
                               unsigned int z) const
  { return 0.25 * (grid_.get_hy(x, y, z) 
                   + grid_.get_hy(x, y + 1, z)
                   + grid_.get_hy(x - 1, y, z)
                   + grid_.get_hy(x - 1, y + 1, z)); }

  inline field_t get_avg_h_t2(unsigned int x, unsigned int y, 
                               unsigned int z) const
  { return 0.25 * (grid_.get_hz(x, y, z) 
                   + grid_.get_hz(x, y, z + 1)
                   + grid_.get_hz(x - 1, y, z)
                   + grid_.get_hz(x - 1, y, z + 1)); }
};


/**
 * Grid adapter class for field components on the Y = constant or XZ
 * plane. 
 *
 * Et1 = EZ
 * Et2 = EX
 * Ht1 = HZ
 * Ht2 = HX
 */
class XZPlane : public GridPlane
{
private:
protected:
public:
  XZPlane(Grid &grid)
    : GridPlane(grid)
  {}

  ~XZPlane()
  {}

  inline const field_t *get_e_n_ptr(unsigned int x, unsigned int y, 
                                    unsigned int z) const
  { return grid_.get_pointer(grid_point(x, y, z), FC_EY); }

  inline const field_t *get_h_n_ptr(unsigned int x, unsigned int y, 
                                    unsigned int z) const
  { return grid_.get_pointer(grid_point(x, y, z), FC_HY); }

  inline const field_t *get_e_t1_ptr(unsigned int x, unsigned int y, 
                                     unsigned int z) const
  { return grid_.get_pointer(grid_point(x, y, z), FC_EZ); }

  inline const field_t *get_h_t1_ptr(unsigned int x, unsigned int y, 
                                     unsigned int z) const
  { return grid_.get_pointer(grid_point(x, y, z), FC_HZ); }

  inline const field_t *get_e_t2_ptr(unsigned int x, unsigned int y, 
                                     unsigned int z) const
  { return grid_.get_pointer(grid_point(x, y, z), FC_EX); }

  inline const field_t *get_h_t2_ptr(unsigned int x, unsigned int y, 
                                     unsigned int z) const
  { return grid_.get_pointer(grid_point(x, y, z), FC_HX); }
  
  inline void set_e_n(unsigned int x, unsigned int y, 
                       unsigned int z, field_t value)
  { grid_.set_ey(x, y, z, value); }

  inline void set_h_n(unsigned int x, unsigned int y, 
                       unsigned int z, field_t value)
  { grid_.set_hy(x, y, z, value); }

  inline void set_e_t1(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value) 
  { grid_.set_ez(x, y, z, value); }

  inline void set_e_t2(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  { grid_.set_ex(x, y, z, value); }

  inline void set_h_t1(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value) 
  { grid_.set_hz(x, y, z, value); }

  inline void set_h_t2(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  { grid_.set_hx(x, y, z, value); }

  inline field_t get_e_n(unsigned int x, unsigned int y, 
                          unsigned int z) const
  { return grid_.get_ey(x, y, z); }

  inline field_t get_h_n(unsigned int x, unsigned int y, 
                          unsigned int z) const
  { return grid_.get_hy(x, y, z); }

  inline field_t get_e_t1(unsigned int x, unsigned int y, 
                           unsigned int z) const
  { return grid_.get_ez(x, y, z); }

  inline field_t get_e_t2(unsigned int x, unsigned int y, 
                           unsigned int z) const
  { return grid_.get_ex(x, y, z); }

  inline field_t get_h_t1(unsigned int x, unsigned int y, 
                           unsigned int z) const
  { return grid_.get_hz(x, y, z); }

  inline field_t get_h_t2(unsigned int x, unsigned int y, 
                           unsigned int z) const
  { return grid_.get_hx(x, y, z); }

  inline field_t get_avg_e_t1(unsigned int x, unsigned int y, 
                               unsigned int z) const
  { return 0.5 * (grid_.get_ez(x, y, z) + grid_.get_ez(x + 1, y, z)); }

  inline field_t get_avg_e_t2(unsigned int x, unsigned int y, 
                               unsigned int z) const
  { return 0.5 * (grid_.get_ex(x, y, z) + grid_.get_ex(x, y, z + 1)); }

  inline field_t get_avg_h_t1(unsigned int x, unsigned int y, 
                               unsigned int z) const
  { return 0.25 * (grid_.get_hz(x, y, z) 
                   + grid_.get_hz(x, y - 1, z)
                   + grid_.get_hz(x, y, z + 1)
                   + grid_.get_hz(x, y - 1, z + 1)); }

  inline field_t get_avg_h_t2(unsigned int x, unsigned int y, 
                               unsigned int z) const
  { return 0.25 * (grid_.get_hx(x, y, z) 
                   + grid_.get_hx(x, y - 1, z)
                   + grid_.get_hx(x + 1, y, z)
                   + grid_.get_hx(x + 1, y - 1, z)); }

};

/**
 * Grid adapter class for field components on the Z = constant or XY
 * plane. 
 *
 * Et1 = EX
 * Et2 = EY
 * Ht1 = HX
 * Ht2 = HY
 */
class XYPlane : public GridPlane
{
private:
protected:
public:
  XYPlane(Grid &grid)
    : GridPlane(grid)
  {}

  ~XYPlane()
  {}
  
  inline const field_t *get_e_n_ptr(unsigned int x, unsigned int y, 
                                    unsigned int z) const
  {
    return grid_.get_pointer(grid_point(x, y, z), FC_EZ);
  }

  inline const field_t *get_h_n_ptr(unsigned int x, unsigned int y, 
                                    unsigned int z) const
  {
    return grid_.get_pointer(grid_point(x, y, z), FC_HZ);
  }

  inline const field_t *get_e_t1_ptr(unsigned int x, unsigned int y, 
                                     unsigned int z) const
  {
    return grid_.get_pointer(grid_point(x, y, z), FC_EX);
  }

  inline const field_t *get_h_t1_ptr(unsigned int x, unsigned int y, 
                                     unsigned int z) const
  {
    return grid_.get_pointer(grid_point(x, y, z), FC_HX);
  }

  inline const field_t *get_e_t2_ptr(unsigned int x, unsigned int y, 
                                     unsigned int z) const
  {
    return grid_.get_pointer(grid_point(x, y, z), FC_EY);
  }

  inline const field_t *get_h_t2_ptr(unsigned int x, unsigned int y, 
                                     unsigned int z) const
  {
    return grid_.get_pointer(grid_point(x, y, z), FC_HY);
  }

  inline void set_e_n(unsigned int x, unsigned int y, 
                       unsigned int z, field_t value)
  {
    grid_.set_ez(x, y, z, value);
  }

  inline void set_h_n(unsigned int x, unsigned int y, 
                       unsigned int z, field_t value)
  {
    grid_.set_hz(x, y, z, value);
  }

  inline void set_e_t1(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  {
    grid_.set_ex(x, y, z, value);
  }

  inline void set_e_t2(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  {
    grid_.set_ey(x, y, z, value);
  }

  inline void set_h_t1(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  {
    grid_.set_hx(x, y, z, value);
  }

  inline void set_h_t2(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  {
    grid_.set_hy(x, y, z, value);
  }

  inline field_t get_e_n(unsigned int x, unsigned int y, 
                          unsigned int z) const
  {
    return grid_.get_ez(x, y, z);
  }

  inline field_t get_h_n(unsigned int x, unsigned int y, 
                          unsigned int z) const
  {
    return grid_.get_hz(x, y, z);
  }

  inline field_t get_e_t1(unsigned int x, unsigned int y, 
                           unsigned int z) const
  {
    return grid_.get_ex(x, y, z);
  }

  inline field_t get_e_t2(unsigned int x, unsigned int y, 
                           unsigned int z) const
  {
    return grid_.get_ey(x, y, z);
  }

  inline field_t get_h_t1(unsigned int x, unsigned int y, 
                           unsigned int z) const
  {
    return grid_.get_hx(x, y, z);
  }

  inline field_t get_h_t2(unsigned int x, unsigned int y, 
                           unsigned int z) const
  {
    return grid_.get_hy(x, y, z);
  }

  inline field_t get_avg_e_t1(unsigned int x, unsigned int y, 
                               unsigned int z) const
  { return 0.5 * (grid_.get_ex(x, y, z) + grid_.get_ex(x, y + 1, z)); }

  inline field_t get_avg_e_t2(unsigned int x, unsigned int y, 
                               unsigned int z) const
  { return 0.5 * (grid_.get_ey(x, y, z) + grid_.get_ey(x + 1, y, z)); }

  inline field_t get_avg_h_t1(unsigned int x, unsigned int y, 
                               unsigned int z) const
  { return 0.25 * (grid_.get_hx(x, y, z) 
                   + grid_.get_hx(x, y, z - 1)
                   + grid_.get_hx(x + 1, y, z)
                   + grid_.get_hx(x + 1, y, z - 1)); }

  inline field_t get_avg_h_t2(unsigned int x, unsigned int y, 
                               unsigned int z) const
  { return 0.25 * (grid_.get_hy(x, y, z) 
                   + grid_.get_hy(x, y, z - 1)
                   + grid_.get_hy(x, y + 1, z)
                   + grid_.get_hy(x, y + 1, z - 1)); }
};

#endif // GRID_FACE_H
