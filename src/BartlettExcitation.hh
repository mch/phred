#ifndef BARTLETT_EXCITATION_H
#define BARTLETT_EXCITATION_H

#include "Excitation.hh"

/**
 * An excitation that is windowed in 3 space by Bartlett windows (tent
 * functions). This should make a pseudo plane wave...
 * 
 * Check the Bartlett Window definition at:
 * http://www.cg.tuwien.ac.at/studentwork/CESCG/CESCG99/TTheussl/node6.html
 */
template<class T>
class BartlettExcitation : public Excitation<T>
{
private:
protected:
public:
  BartlettExcitation(T &sf) 
    : Excitation<T>(sf)
  {}

  ~BartlettExcitation()
  {}

  /**
   * This applied the windowed excitation to the grid.
   *
   * @param grid the grid to mess with
   * @param time_step the time step we are on
   * @param type the field components to excite (E or H)
   */
  void excite(Grid &grid, unsigned int time_step, 
              FieldType type)
  {
    // Find out where we fit in this grid (convert to local coordinates)
    if (type != BOTH && type != type_)
      return;

    region_t r = grid.global_to_local(region_);

    field_t sf = sf_.source_function(grid, time_step);
    field_t fld[3];

    // Bartlett window coefficients
    field_t wx = 0, wy = 0, wz = 0, w = 1;

    fld[0] = sf * polarization_[0];
    fld[1] = sf * polarization_[1];
    fld[2] = sf * polarization_[2];

    if (!soft_) 
    {
      for(unsigned int i = r.xmin; i < r.xmax; i++)
      {
        wx = 1 - fabs( ( (i - r.xmin) - 0.5 * (r.xmax - r.xmin - 1)) 
                       / (0.5 * (r.xmax - r.xmin + 1)));

        for (unsigned int j = r.ymin; j < r.ymax; j++)
        {
          wy = 1 - fabs( ( (j - r.ymin) - 0.5 * (r.ymax - r.ymin - 1)) 
                         / (0.5 * (r.ymax - r.ymin + 1)));

          for (unsigned int k = r.zmin; k < r.zmax; k++)
          {
            wz = 1 - fabs( ( (k - r.zmin) - 0.5 * (r.zmax - r.zmin - 1)) 
                           / (0.5 * (r.zmax - r.zmin + 1)));

            w = wx * wy * wz;

            switch (type_) 
            {
            case E:
              if (polarization_[0] != 0.0) grid.set_ex(i,j,k, fld[0] * w);
              if (polarization_[1] != 0.0) grid.set_ey(i,j,k, fld[1] * w);
              if (polarization_[2] != 0.0) grid.set_ez(i,j,k, fld[2] * w);
              break;

            case H:
              if (polarization_[0] != 0.0) grid.set_hx(i,j,k, fld[0] * w);
              if (polarization_[1] != 0.0) grid.set_hy(i,j,k, fld[1] * w);
              if (polarization_[2] != 0.0) grid.set_hz(i,j,k, fld[2] * w);
              break;

            case BOTH: // Isn't meant for Excitations.
              throw std::exception();
              break; 
            }
          }
        }
      }
    } else {
      for(unsigned int i = r.xmin; i < r.xmax; i++)
      {
        wx = 1 - fabs( ( (i - r.xmin) - 0.5 * (r.xmax - r.xmin - 1)) 
                       / (0.5 * (r.xmax - r.xmin + 1)));

        for (unsigned int j = r.ymin; j < r.ymax; j++)
        {
          wy = 1 - fabs( ( (j - r.ymin) - 0.5 * (r.ymax - r.ymin - 1)) 
                         / (0.5 * (r.ymax - r.ymin + 1)));

          for (unsigned int k = r.zmin; k < r.zmax; k++)
          {
            wz = 1 - fabs( ( (k - r.zmin) - 0.5 * (r.zmax - r.zmin - 1)) 
                           / (0.5 * (r.zmax - r.zmin + 1)));

            w = wx * wy * wz;

            switch (type_) 
            {
            case E:
              if (polarization_[0] != 0.0) grid.set_ex(i,j,k, 
                                                       w * fld[0] 
                                                       + grid.get_ex(i,j,k));
              if (polarization_[1] != 0.0) grid.set_ey(i,j,k, 
                                                       w * fld[1] 
                                                       + grid.get_ey(i,j,k));
              if (polarization_[2] != 0.0) grid.set_ez(i,j,k, 
                                                       w * fld[2] 
                                                       + grid.get_ez(i,j,k));
              break;

            case H:
              if (polarization_[0] != 0.0) grid.set_hx(i,j,k, 
                                                       w * fld[0] 
                                                       + grid.get_hx(i,j,k));
              if (polarization_[1] != 0.0) grid.set_hy(i,j,k, 
                                                       w * fld[1] 
                                                       + grid.get_hy(i,j,k));
              if (polarization_[2] != 0.0) grid.set_hz(i,j,k, 
                                                       w * fld[2] 
                                                       + grid.get_hz(i,j,k));
              break;

            case BOTH: // Isn't meant for Excitations.
              throw std::exception();
              break; 
            }
          }
        }
      }
    }
  }
  
};

#endif // BARTLETT_EXCITATION_H
