#ifndef GEOMETRY_H
#define GEOMETRY_H

/**
 * This is an abstract base class which defines an interface which is
 * to be used to build objects which define... geometry
 * objects. Boxes and spheres and stuff. You know. Each object must
 * know how to iterate over the cells it constains, so that you can
 * do things to it. 
 */
class Geometry
{
private:
protected:
public:
  Geometry();
  virtual ~Geometry() = 0;

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
    const iterator &operator++(const iterator &rhs);
    
    /**
     * Dereference the iterator; return a index into the grid where
     * the current element can be found. 
     */
    unsigned int operator*(const iterator &rhs);
  };
};

#endif // GEOMETRY_H
