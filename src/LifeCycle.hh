#ifndef LIFE_CYCLE_H
#define LIFE_CYCLE_H

#include "Grid.hh"

/**
 * An interface that defines life cycle functions used by a FDTD sim
 * management object to set up member objects when leaving define mode
 * so that they can allocate memory and init constants that may depend
 * on the state of the grid or other objects. 
 */
class LifeCycle
{
private:
protected:
public:
  LifeCycle();

  virtual ~LifeCycle() = 0;

  /**
   * Subclasses can implement this to allocate memory or init
   * constants or whatever. Called just before the simulation
   * starts. The default implementation does nothing. 
   */
  void init(const Grid &grid)
  {}
  
  /**
   * Subclasses can implement this to deallocate memory or
   * whatever. Called just after the simulation ends. The default
   * implementation does nothing. 
   */
  void deinit(const Grid &grid)
  {}
};

#endif // LIFE_CYCLE_H
