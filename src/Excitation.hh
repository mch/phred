// An abstract base class for excitations. Subclasses contain the data
// and an implementation of excite(Grid &grid) which knows how to
// apply the excitation to the grid. 

#ifndef EXCITE_H
#define EXCITE_H

#include "Types.hh"
#include "Grid.hh"

class Excitation 
{
protected:
  // Region to apply the source to (in global coordinates)
  unsigned int x_start_;
  unsigned int y_start_;
  unsigned int z_start_;

  unsigned int x_end_;
  unsigned int y_end_;
  unsigned int z_end_;

  // Field component to excite
  FieldComponent component_;

public:
  Excitation();
  virtual ~Excitation();

  /**
   * This function is called to apply this excitation to the grid. 
   *
   * @param the grid to operate on
   * @param the time step we are on
   */
  virtual void excite(Grid &grid, unsigned int time_step) = 0;

  /**
   * Set the (rectangular) region to apply the excitation to
   *
   * @param The starting x coordinate
   * @param The starting y coordinate
   * @param The starting z coordinate
   *
   * @param The ending x coordinate
   * @param The ending y coordinate
   * @param The ending z coordinate
   */ 
  void set_region(unsigned int x_start, unsigned int y_start, 
                  unsigned int z_start, unsigned int x_end, 
                  unsigned int y_end, unsigned int z_end);
                  
};

#endif // EXCITE_H

