#ifndef BOX_H_
#define BOX_H_

#include "Geometry.hh"

/**
 * Defines a box in the grid
 */
class Box : public Geometry
{
private:
protected:
  region_t r_;

public:
  Box();
  Box(region_t r);
  ~Box();

  /**
   * Set the region to operate on
   */
  inline void set_region(region_t r)
  {
    r_ = r;
  }

  /**
   * Set the region to operate on
   */
  void set_region(unsigned int xstart, unsigned int xstop, 
                  unsigned int ystart, unsigned int ystop, 
                  unsigned int zstart, unsigned int zstop);

  /**
   * Init the region, set the material in the grid. Normally called by
   * the FDTD object.
   */
  virtual void init(const Grid &grid);
  
  /**
   * Update the material indicies of the grid 
   */
  virtual void set_material(Grid &grid);

};

#endif
