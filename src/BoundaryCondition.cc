#include "BoundaryCondition.hh"

BoundaryCond::region BoundaryCond::find_face(Face face, Grid &grid)
{
  region r;

  switch (face) 
  {
  case FRONT:
    r.xmax = r.xmin = grid.get_ldx();
    r.ymin = r.zmin = 0;
    r.ymax = grid.get_ldy();
    r.zmax = grid.get_ldz();
    break;

  case BACK:
    r.xmax = r.xmin = 0;
    r.ymin = r.zmin = 0;
    r.ymax = grid.get_ldy();
    r.zmax = grid.get_ldz();
    break;

  case TOP:
    r.zmin = r.zmax = grid.get_ldz();
    r.ymin = r.xmin = 0;
    r.xmax = grid.get_ldx();
    r.ymax = grid.get_ldy();
    break;

  case BOTTOM:
    r.zmin = r.zmax = 0;
    r.ymin = r.xmin = 0;
    r.xmax = grid.get_ldx();
    r.ymax = grid.get_ldy();
    break;

  case LEFT:
    r.ymin = r.ymax = grid.get_ldy();
    r.zmin = r.xmin = 0;
    r.xmax = grid.get_ldx();
    r.zmax = grid.get_ldz();
    break;

  case RIGHT:
    r.ymin = r.ymax = grid.get_ldy();
    r.zmin = r.xmin = 0;
    r.xmax = grid.get_ldx();
    r.zmax = grid.get_ldz();
    break;
  }

  return r; 
}
