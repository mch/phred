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
 * component.
 */
class GridPlane
{
private:
protected:
  Grid &grid_;

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
                          unsigned int z) = 0;

  /** 
   * Get the normal H field component
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  virtual field_t get_h_n(unsigned int x, unsigned int y, 
                          unsigned int z) = 0;

  /**
   * Get the first tangential E field component, E_1. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  virtual field_t get_e_t1(unsigned int x, unsigned int y, 
                           unsigned int z) = 0;

  /**
   * Get the second tangential E field component, E_2. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  virtual field_t get_e_t2(unsigned int x, unsigned int y, 
                           unsigned int z) = 0;

  /**
   * Get the first tangential H field component, H_1. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  virtual field_t get_h_t1(unsigned int x, unsigned int y, 
                           unsigned int z) = 0;

  /**
   * Get the second tangential H field component, H_2. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  virtual field_t get_h_t2(unsigned int x, unsigned int y, 
                           unsigned int z) = 0;

};

/**
 * Grid adapter class for field components on the X = constant or YZ
 * plane. 
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
  
  /** 
   * Set the normal E field component, Ex
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  inline void set_e_n(unsigned int x, unsigned int y, 
                       unsigned int z, field_t value)
  {
    grid_.set_ex(x, y, z, value);
  }

  /** 
   * Set the normal H field component, Hx
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  inline void set_h_n(unsigned int x, unsigned int y, 
                       unsigned int z, field_t value)
  {
    grid_.set_hx(x, y, z, value);
  }

  /**
   * Set the first tangential E field component, Ey.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  inline void set_e_t1(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  {
    grid_.set_ey(x, y, z, value);
  }

  /**
   * Set the second tangential E field component, Ez. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  inline void set_e_t2(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  {
    grid_.set_ez(x, y, z, value);
  }

  /**
   * Set the first tangential H field component, Hy.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  inline void set_h_t1(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value) 
  {
    grid_.set_hy(x, y, z, value);
  }

  /**
   * Set the second tangential H field component, Hz.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  inline void set_h_t2(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  {
    grid_.set_hz(x, y, z, value);
  }

  /** 
   * Get the normal E field component, Ex
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  inline field_t get_e_n(unsigned int x, unsigned int y, 
                          unsigned int z)
  {
    return grid_.get_ex(x, y, z);
  }

  /** 
   * Get the normal H field component, Hx
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  inline field_t get_h_n(unsigned int x, unsigned int y, 
                          unsigned int z)
  {
    return grid_.get_hx(x, y, z);
  }

  /**
   * Get the first tangential E field component, Ey.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  inline field_t get_e_t1(unsigned int x, unsigned int y, 
                           unsigned int z)
  {
    return grid_.get_ey(x, y, z);
  }

  /**
   * Get the second tangential E field component, Ez.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  inline field_t get_e_t2(unsigned int x, unsigned int y, 
                           unsigned int z)
  {
    return grid_.get_ez(x, y, z);
  }

  /**
   * Get the first tangential H field component, Hy.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  inline field_t get_h_t1(unsigned int x, unsigned int y, 
                           unsigned int z)
  {
    return grid_.get_hy(x, y, z);
  }

  /**
   * Get the second tangential H field component, Hz.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  inline field_t get_h_t2(unsigned int x, unsigned int y, 
                           unsigned int z)
  {
    return grid_.get_hz(x, y, z);
  }
};


/**
 * Grid adapter class for field components on the Y = constant or XZ
 * plane. 
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
  
  /** 
   * Set the normal E field component, Ey
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  inline void set_e_n(unsigned int x, unsigned int y, 
                       unsigned int z, field_t value)
  {
    grid_.set_ey(x, y, z, value);
  }

  /** 
   * Set the normal H field component, Hy
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  inline void set_h_n(unsigned int x, unsigned int y, 
                       unsigned int z, field_t value)
  {
    grid_.set_hy(x, y, z, value);
  }

  /**
   * Set the first tangential E field component, Ez.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  inline void set_e_t1(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  {
    grid_.set_ez(x, y, z, value);
  }

  /**
   * Set the second tangential E field component, Ex. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  inline void set_e_t2(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  {
    grid_.set_ex(x, y, z, value);
  }

  /**
   * Set the first tangential H field component, Hz.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  inline void set_h_t1(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value) 
  {
    grid_.set_hz(x, y, z, value);
  }

  /**
   * Set the second tangential H field component, Hx.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  inline void set_h_t2(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  {
    grid_.set_hx(x, y, z, value);
  }

  /** 
   * Get the normal E field component, Ey
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  inline field_t get_e_n(unsigned int x, unsigned int y, 
                          unsigned int z)
  {
    return grid_.get_ey(x, y, z);
  }

  /** 
   * Get the normal H field component, Hy
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  inline field_t get_h_n(unsigned int x, unsigned int y, 
                          unsigned int z)
  {
    return grid_.get_hy(x, y, z);
  }

  /**
   * Get the first tangential E field component, Ez.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  inline field_t get_e_t1(unsigned int x, unsigned int y, 
                           unsigned int z)
  {
    return grid_.get_ez(x, y, z);
  }

  /**
   * Get the second tangential E field component, Ex.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  inline field_t get_e_t2(unsigned int x, unsigned int y, 
                           unsigned int z)
  {
    return grid_.get_ex(x, y, z);
  }

  /**
   * Get the first tangential H field component, Hz.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  inline field_t get_h_t1(unsigned int x, unsigned int y, 
                           unsigned int z)
  {
    return grid_.get_hz(x, y, z);
  }

  /**
   * Get the second tangential H field component, Hx.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  inline field_t get_h_t2(unsigned int x, unsigned int y, 
                           unsigned int z)
  {
    return grid_.get_hx(x, y, z);
  }
};

/**
 * Grid adapter class for field components on the Z = constant or XY
 * plane. 
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
  
  /** 
   * Set the normal E field component, Ez
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  inline void set_e_n(unsigned int x, unsigned int y, 
                       unsigned int z, field_t value)
  {
    grid_.set_ez(x, y, z, value);
  }

  /** 
   * Set the normal H field component, Hz
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  inline void set_h_n(unsigned int x, unsigned int y, 
                       unsigned int z, field_t value)
  {
    grid_.set_hz(x, y, z, value);
  }

  /**
   * Set the first tangential E field component, Ex.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  inline void set_e_t1(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  {
    grid_.set_ex(x, y, z, value);
  }

  /**
   * Set the second tangential E field component, Ey. 
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  inline void set_e_t2(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  {
    grid_.set_ey(x, y, z, value);
  }

  /**
   * Set the first tangential H field component, Hx.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  inline void set_h_t1(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value) 
  {
    grid_.set_hx(x, y, z, value);
  }

  /**
   * Set the second tangential H field component, Hy.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param value new field component value
   */
  inline void set_h_t2(unsigned int x, unsigned int y, 
                        unsigned int z, field_t value)
  {
    grid_.set_hy(x, y, z, value);
  }

  /** 
   * Get the normal E field component, Ez
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  inline field_t get_e_n(unsigned int x, unsigned int y, 
                          unsigned int z)
  {
    return grid_.get_ez(x, y, z);
  }

  /** 
   * Get the normal H field component, Hz
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  inline field_t get_h_n(unsigned int x, unsigned int y, 
                          unsigned int z)
  {
    return grid_.get_hz(x, y, z);
  }

  /**
   * Get the first tangential E field component, Ex.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  inline field_t get_e_t1(unsigned int x, unsigned int y, 
                           unsigned int z)
  {
    return grid_.get_ex(x, y, z);
  }

  /**
   * Get the second tangential E field component, Ey.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  inline field_t get_e_t2(unsigned int x, unsigned int y, 
                           unsigned int z)
  {
    return grid_.get_ey(x, y, z);
  }

  /**
   * Get the first tangential H field component, Hx.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  inline field_t get_h_t1(unsigned int x, unsigned int y, 
                           unsigned int z)
  {
    return grid_.get_hx(x, y, z);
  }

  /**
   * Get the second tangential H field component, Hy.
   *
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @return field component value
   */
  inline field_t get_h_t2(unsigned int x, unsigned int y, 
                           unsigned int z)
  {
    return grid_.get_hy(x, y, z);
  }
};

#endif // GRID_FACE_H
