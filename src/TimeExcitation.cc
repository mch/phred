#include "TimeExcitation.hh"
#include <exception>

void TimeExcitation::excite(Grid &grid, unsigned int time_step,
                            FieldType type)
{
  // Find out where we fit in this grid (convert to local coordinates)
  if (type != BOTH && type == type_)
    return;

  region_t r = grid.global_to_local(x_start_, x_end_, 
                                    y_start_, y_end_,
                                    z_start_, z_end_);

  field_t sf = source_function(grid, time_step);
  field_t fld[3];

  fld[0] = sf * polarization_[0];
  fld[1] = sf * polarization_[1];
  fld[2] = sf * polarization_[2];

  for(unsigned int i = r.xmin; i < r.xmax; i++)
  {
    for (unsigned int j = r.ymin; j < r.ymax; j++)
    {
      for (unsigned int k = r.zmin; k < r.zmax; k++)
      {
        switch (type_) 
        {
        case E:
          if (fld[0] > 0) grid.set_ex(i,j,k, fld[0]);
          if (fld[1] > 0) grid.set_ey(i,j,k, fld[1]);
          if (fld[2] > 0) grid.set_ez(i,j,k, fld[2]);
          break;

        case H:
          if (fld[0] > 0) grid.set_hx(i,j,k, fld[0]);
          if (fld[1] > 0) grid.set_hy(i,j,k, fld[1]);
          if (fld[2] > 0) grid.set_hz(i,j,k, fld[2]);
          break;

        case BOTH: // Isn't meant for Excitations.
          throw std::exception();
          break; 
        }
      }
    }
  }
}
