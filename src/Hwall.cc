#include "Hwall.hh"
#include "GridPlane.hh"

Hwall::Hwall()
{}

Hwall::~Hwall()
{}

template <class T>
void Hwall::condition(region r, Grid &grid)
{
  T p(grid);

  for (unsigned int i = r.xmin; i < r.xmax; i++) {
    for (unsigned int j = r.ymin; j < r.ymax; j++) {
      for (unsigned int k = r.zmin; k < r.zmax; k++) {
        p.set_h_t1(i, j, k, 0.);
        p.set_h_t2(i, j, k, 0.);
      }
    }
  }
}

void Hwall::apply(Face face, Grid &grid)
{
  region r = find_face(face, grid);

  switch (face)
  {
  case FRONT:
  case BACK:
    condition<YZPlane>(r, grid);
    break;

  case LEFT:
  case RIGHT:
    condition<XZPlane>(r, grid);
    break;

  case TOP:
  case BOTTOM:
    condition<XYPlane>(r, grid);
    break;
  }
  
}




