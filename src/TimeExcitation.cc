#include "TimeExcitation.hh"

void TimeExcitation::excite(Grid &grid, unsigned int time_step)
{
  // Find out where we fit in this grid (convert to local coordinates)

  region r = grid.global_to_local(x_start_, x_end_, 
                                  y_start_, y_end_,
                                  z_start_, z_end_);

  field_t sf = source_function(grid, time_step);

  for(unsigned int i = r.xmin; i < r.xmax; i++)
  {
    for (unsigned int j = r.ymin; j < r.ymax; j++)
    {
      for (unsigned int k = r.zmin; k < r.zmax; k++)
      {
        switch (type_) 
        {
        case E:
          grid.set_ex(i,j,k, sf * polarization_[0]);
          grid.set_ey(i,j,k, sf * polarization_[1]);
          grid.set_ez(i,j,k, sf * polarization_[2]);
          break;

        case H:
          grid.set_hx(i,j,k, sf * polarization_[0]);
          grid.set_hy(i,j,k, sf * polarization_[1]);
          grid.set_hz(i,j,k, sf * polarization_[2]);
          break;
        }
      }
    }
  }
}
