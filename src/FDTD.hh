#ifndef FDTD_H
#define FDTD_H

#include <map>
#include <vector>
#include <string>

using namespace std;

#include "Types.hh"

/**
 * This is a sort of convience wrapper object that runs the
 * simulation. It is primarily intended to be used from python. It
 * holds lists of results, excitations, data writers, etc, manages
 * thier memory and life cycles, etc. 
 */
class FDTD
{
private:
protected:
public:
  FDTD();
  virtual ~FDTD();

  /**
   * Set the global (entire problem) grid size
   */
  void set_grid_size(unsigned int x, 
                     unsigned int y, unsigned int z);

  /**
   * Set the grid deltas, the size of the cells. The time delta is
   * automatically computed using the stability condition. 
   */
  void set_grid_deltas(field_t dx, field_t dy, field_t dz);

  /**
   * 
};

#endif // FDTD_H
