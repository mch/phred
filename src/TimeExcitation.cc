#include "TimeExcitation.hh"

void TimeExcitation::excite(Grid &grid, unsigned int time_step)
{
  // Find out where we fit in this grid (convert to local coordinates)
  unsigned int xs, ys, zs, xe, ye, ze;

  xs = (grid.get_lsx() > x_start_) ? grid.get_lsx() 
    : x_start_ - grid.get_lsx();
  ys = (grid.get_lsy() > y_start_) ? grid.get_lsy() 
    : y_start_ - grid.get_lsy();
  zs = (grid.get_lsz() > z_start_) ? grid.get_lsz() 
    : z_start_ - grid.get_lsz();

  xe = (grid.get_lsx() + grid.get_gdx() > x_end_) 
    ? grid.get_lsx() + grid.get_gdx() : x_end_ - grid.get_lsx();
  ye = (grid.get_lsy() + grid.get_gdy() > y_end_) 
    ? grid.get_lsy() + grid.get_gdy() : y_end_ - grid.get_lsy();
  ze = (grid.get_lsz() + grid.get_gdz() > z_end_) 
    ? grid.get_lsz() + grid.get_gdz() : z_end_ - grid.get_lsz();

  field_t sf = source_function(grid, time_step);

  for(unsigned int i = xs; i < xe; i++)
  {
    for (unsigned int j = ys; j < ye; j++)
    {
      for (unsigned int k = zs; k < ze; k++)
      {
        switch (component_) 
        {
        case EX:
          grid.set_ex(i,j,k, sf);
          break;

        case EY:
          grid.set_ey(i,j,k, sf);          
          break;

        case EZ:
          grid.set_ez(i,j,k, sf);
          break;

        case HX:
          grid.set_hx(i,j,k, sf);          
          break;

        case HY:
          grid.set_hy(i,j,k, sf);                    
          break;

        case HZ:
          grid.set_hz(i,j,k, sf);          
          break;
        }
      }
    }
  }
}
