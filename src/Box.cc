#include "Box.hh"
#include "Grid.hh"

Box::Box()
{}

Box::~Box()
{}

void Box::set_region(unsigned int xstart, unsigned int xstop, 
                     unsigned int ystart, unsigned int ystop, 
                     unsigned int zstart, unsigned int zstop)
{
  r_.xmin = xstart;
  r_.xmax = xstop;
  r_.ymin = ystart;
  r_.ymax = ystop;
  r_.zmin = zstart;
  r_.zmax = zstop;
}

void Box::init(const Grid &grid)
{

}
void Box::set_material(Grid &grid)
{
  region_t r = grid.global_to_local(r_);

  for (unsigned int i = r.xmin; i < r.xmax; i++)
  {
    for (unsigned int j = r.ymin; j < r.ymax; j++)
    {
      for (unsigned int k = r.zmin; k < r.zmax; k++)
      {
        grid.set_material(i, j, k, material_id_);
      }
    }
  }
}
