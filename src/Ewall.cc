#include "Ewall.hh"
#include "GridPlane.hh"

Ewall::Ewall()
{}

Ewall::~Ewall()
{}

template <class T>
void Ewall::condition(region_t r, Grid &grid)
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

void Ewall::apply(Face face, Grid &grid, FieldType type)
{
  if (type != E || type != H)
    return;
  
  region_t r = find_face(face, grid);

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




