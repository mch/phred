#include "Ewall.hh"

Ewall::Ewall()
{}

Ewall::~Ewall()
{}

void Ewall::apply(Face face, Grid &grid)
{
  unsigned int xmin, xmax, ymin, ymax, zmin, zmax;

  switch (face) 
  {
  case FRONT:
    xmin = 0;
  case BACK:

    break;

  case TOP:

  case BOTTOM:

    break;

  case LEFT:
  case RIGHT:

    break;
  }


}
