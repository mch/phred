#include "Ewall.hh"

Ewall::Ewall()
{}

Ewall::~Ewall()
{}

void Ewall::apply(Face face, Grid &grid)
{
  region r = find_face(face, grid);

  switch (face)
  {
  case FRONT:
  case BACK:
    helper<YZPlane>(r, grid);
    break;

  case LEFT:
  case RIGHT:
    helper<XZPlane>(r, grid);
    break;

  case TOP:
  case BOTTOM:
    helper<XYPlane>(r, grid);
    break;
  }
  
}

template <class T>
void Ewall::helper(region r, Grid &grid)
{
  T p(grid);

  for (unsigned int i = r.xmin; i < r.xmax; i++) {
    for (unsigned int j = r.ymin; j < r.ymax; j++) {
      for (unsigned int k = r.zmin; k < r.zmax; k++) {
        p.set_e_t1(i, j, k, 0.);
        p.set_e_t2(i, j, k, 0.);
      }
    }
  }  
}
