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


};

#endif // GEOMETRY_H
