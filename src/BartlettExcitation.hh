#ifndef BARTLETT_EXCITATION_H
#define BARTLETT_EXCITATION_H

#include "WindowedExcitation.hh"

/**
 * An excitation that is windowed in 3 space by Bartlett windows (tent
 * functions). This should make a pseudo plane wave...
 * 
 * Check the Bartlett Window definition at:
 * http://www.cg.tuwien.ac.at/studentwork/CESCG/CESCG99/TTheussl/node6.html
 */
class BartlettExcitation : public WindowedExcitation
{
private:
protected:
  /**
   * Defines a bartlett window function. 
   */
  virtual field_t window(region_t r, unsigned int i, unsigned int j, 
                         unsigned int k) 
  {
    field_t wx = 0, wy = 0, wz = 0;

    wx = 1 - fabs( ( (i - r.xmin) - 0.5 * (r.xmax - r.xmin - 1)) 
                   / (0.5 * (r.xmax - r.xmin + 1)));
    
    wy = 1 - fabs( ( (j - r.ymin) - 0.5 * (r.ymax - r.ymin - 1)) 
                   / (0.5 * (r.ymax - r.ymin + 1)));

    wz = 1 - fabs( ( (k - r.zmin) - 0.5 * (r.zmax - r.zmin - 1)) 
                   / (0.5 * (r.zmax - r.zmin + 1)));
    
    return wx * wy * wz;
  }

public:
  BartlettExcitation(SourceFunction *sf) 
    : WindowedExcitation(sf)
  {}

  ~BartlettExcitation()
  {}
  
};

#endif // BARTLETT_EXCITATION_H
