#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "LifeCycle.hh"
#include "Types.hh"

/**
 * This is an abstract base class which defines an interface which is
 * to be used to build objects which define... geometry
 * objects. Boxes and spheres and stuff. You know. Each object must
 * know how to iterate over the cells it constains, so that you can
 * do things to it. 
 */
class Geometry : public LifeCycle
{
private:
protected:
  unsigned int material_id_;

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
   * inside the grid and fill in the material values of the grid. 
   */ 
  virtual void init(const Grid &grid) = 0;
};

#endif // GEOMETRY_H
