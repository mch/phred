#include "BoundaryCondition.hh"
#include "Grid.hh" // Flesh out the forward declaration

region BoundaryCond::find_face(Face face, Grid &grid)
{
  region r;

  switch (face) 
  {
  case FRONT:
    r.xmax = grid.get_ldx();
    r.xmin = r.xmax - 1;
    r.ymin = r.zmin = 0;
    r.ymax = grid.get_ldy();
    r.zmax = grid.get_ldz();
    break;

  case BACK:
    r.xmax = 1;
    r.xmin = 0;
    r.ymin = r.zmin = 0;
    r.ymax = grid.get_ldy();
    r.zmax = grid.get_ldz();
    break;

  case TOP:
    r.zmin = grid.get_ldz() - 1;
    r.zmax = grid.get_ldz();
    r.ymin = r.xmin = 0;
    r.xmax = grid.get_ldx();
    r.ymax = grid.get_ldy();
    break;

  case BOTTOM:
    r.zmin = 0;
    r.zmax = 1;
    r.ymin = r.xmin = 0;
    r.xmax = grid.get_ldx();
    r.ymax = grid.get_ldy();
    break;

  case LEFT:
    r.ymin = 0;
    r.ymax = 1;
    r.zmin = r.xmin = 0;
    r.xmax = grid.get_ldx();
    r.zmax = grid.get_ldz();
    break;

  case RIGHT:
    r.ymin = grid.get_ldy() - 1;
    r.ymax = grid.get_ldy();
    r.zmin = r.xmin = 0;
    r.xmax = grid.get_ldx();
    r.zmax = grid.get_ldz();
    break;
  }

  return r; 
}

// void BoundaryCond::apply(Face face, Grid &grid)
// {
//   region r = find_face(face, grid);

//   switch (face)
//   {
//   case FRONT:
//   case BACK:
//     condition<YZPlane>(r, grid);
//     break;

//   case LEFT:
//   case RIGHT:
//     condition<XZPlane>(r, grid);
//     break;

//   case TOP:
//   case BOTTOM:
//     helper<XYPlane>(r, grid);
//     break;
//   }
  
// }

unsigned int BoundaryCond::get_thickness()
{
  return thickness_;
}
