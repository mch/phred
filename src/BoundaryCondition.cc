#include "BoundaryCondition.hh"
#include "Grid.hh" // Flesh out the forward declaration

// ARGH, seems like the PML is overlapping the Grid update equations,
// and wreaking havoc. This function seems to return things that are
// one item too thick.

region_t BoundaryCond::find_face(Face face, Grid &grid)
{
  region_t r;
  unsigned int temp = thickness_;
  
  //if (temp > 0)
  //  temp--;

  switch (face) 
  {
  case FRONT:
    r.xmax = grid.get_ldx();
    r.xmin = r.xmax - 1 - temp;
    r.ymin = r.zmin = 0;
    r.ymax = grid.get_ldy();
    r.zmax = grid.get_ldz();
    break;

  case BACK:
    r.xmax = 1 + temp;
    r.xmin = 0;
    r.ymin = r.zmin = 0;
    r.ymax = grid.get_ldy();
    r.zmax = grid.get_ldz();
    break;

  case TOP:
    r.zmin = grid.get_ldz() - 1 - temp;
    r.zmax = grid.get_ldz();
    r.ymin = r.xmin = 0;
    r.xmax = grid.get_ldx();
    r.ymax = grid.get_ldy();
    break;

  case BOTTOM:
    r.zmin = 0;
    r.zmax = 1 + temp;
    r.ymin = r.xmin = 0;
    r.xmax = grid.get_ldx();
    r.ymax = grid.get_ldy();
    break;

  case LEFT:
    r.ymin = 0;
    r.ymax = 1 + temp;
    r.zmin = r.xmin = 0;
    r.xmax = grid.get_ldx();
    r.zmax = grid.get_ldz();
    break;

  case RIGHT:
    r.ymin = grid.get_ldy() - 1 - temp;
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
//   region_t r = find_face(face, grid);

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
