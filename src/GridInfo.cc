#include "GridInfo.hh"

GridInfo::GridInfo() 
  : global_dimx_(0), global_dimy_(0), global_dimz_(0), 
    dimx_(0), dimy_(0), dimz_(0), 
    deltax_(0), deltay_(0), deltaz_(0), deltat_(0)
{
  for (int i = 0; i < 6; i++) {
    face_bc_[i] = 0; 
    //face_ptr_owner_[i] = false;
  }
}

GridInfo::GridInfo(const GridInfo &info) {
  *this = info;
}

GridInfo::~GridInfo()
{
  for (int i = 0; i < 6; i++)
  {
    //if (face_ptr_owner_[i]) {
      delete face_bc_[i];
      face_bc_[i] = 0;
      //}
  }
}

GridInfo& GridInfo::operator=(const GridInfo &info)
{
  if (this == &info) return *this;

  global_dimx_ = info.global_dimx_;
  global_dimy_ = info.global_dimy_;
  global_dimz_ = info.global_dimz_;
  
  start_x_ = info.start_x_;
  start_y_ = info.start_y_;
  start_z_ = info.start_z_;
  
  dimx_ = info.dimx_;
  dimy_ = info.dimy_;
  dimz_ = info.dimz_;
  
  deltax_ = info.deltax_;
  deltay_ = info.deltay_;
  deltaz_ = info.deltaz_;
  deltat_ = info.deltat_;
  
  for (int i = 0; i < 6; i++) {
    face_bc_[i] = info.face_bc_[i];
    //face_ptr_owner_ = false;
  }
  
  return *this;
}

void GridInfo::set_boundary(Face face, BoundaryCond *bc, bool take_ownership)
{
  if (take_ownership)
    //face_ptr_owner_[face] = true;
    face_bc_[face] = new counted_ptr<BoundaryCond>(bc);
  else
    face_bc_[face] = new counted_ptr<BoundaryCond>(bc, 2);
}

// BoundaryCond *GridInfo::copy_bc(BoundaryCond *bc, BoundaryCondition bc_type)
// {
//   BoundaryCond *ret = 0;

//   switch (bc_type) {
//   case SUBDOMAIN:
//     ret = new SubdomainBc(dynamic_cast<const SubdomainBc&>(*bc));
//     break;

//   case EWALL:
//     ret = new Ewall(dynamic_cast<const Ewall&>(*bc));
//     break;

//   case HWALL:
//     ret = new Hwall(dynamic_cast<const Hwall&>(*bc));
//     break;

//   case PML:
//     ret = new Pml(dynamic_cast<const Pml&>(*bc));
//     break;

//   case UNKNOWN:
//     ret = new UnknownBc(dynamic_cast<const UnknownBc&>(*bc));
//     break;
//   }

//   return ret;
// }

unsigned int GridInfo::get_face_thickness(Face face)
{
  unsigned int ret = 0;

  if (face_bc_[face]) {
    ret = face_bc_[face]->get()->get_thickness();
  }

  return ret;
}

void GridInfo::apply_boundaries(Grid &grid, FieldType type)
{
  for (unsigned int i = 0; i < 6; i++)
  {
    if (face_bc_[i]) {
      face_bc_[i]->get()->apply(static_cast<Face>(i), grid, type);
    }
  }
}

