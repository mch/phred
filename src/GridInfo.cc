#include "GridInfo.hh"

GridInfo::GridInfo() 
  : global_dimx_(0), global_dimy_(0), global_dimz_(0), 
    dimx_(0), dimy_(0), dimz_(0), 
    deltax_(0), deltay_(0), deltaz_(0), deltat_(0)
{
  for (int i = 0; i < 6; i++) {
    face_bc_[i] = new UnknownBc();
    face_bc_type_[i] = UNKNOWN;
  }
}

GridInfo::GridInfo(const GridInfo &info) {
  *this = info;
}

GridInfo::~GridInfo()
{
  for (int i = 0; i < 6; i++)
  {
    if (face_bc_[i]) {
      delete face_bc_[i];
      face_bc_[i] = 0;
    }
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
    if (info.face_bc_[i]) {
      face_bc_[i] = copy_bc(info.face_bc_[i], info.face_bc_type_[i]);
      face_bc_type_[i] = info.face_bc_type_[i];
    } else {
      face_bc_[i] = new UnknownBc();
      face_bc_type_[i] = UNKNOWN;
    }
  }
  
  return *this;
}

Pml *GridInfo::set_pml_boundary(Face face, unsigned int thickness, 
                                PmlVariation_t var, float nrml_refl)
{
  Pml *pml = new Pml(var, nrml_refl);
  pml->set_thickness(thickness);
  face_bc_[face] = pml;
  face_bc_type_[face] = PML;

  return pml;
}

BoundaryCond *GridInfo::set_boundary(Face face,
                                     BoundaryCondition bc)
{
  face_bc_type_[face] = bc;

  if (face_bc_[face]) {
    delete face_bc_[face];
  }

  switch (bc) {
  case SUBDOMAIN:
    face_bc_[face] = new SubdomainBc();
    break;

  case EWALL:
    face_bc_[face] = new Ewall();
    break;

  case HWALL:
    face_bc_[face] = new Hwall();
    break;

  case PML:
    face_bc_[face] = new Pml();
    break;

  default:
    face_bc_[face] = new UnknownBc();
    face_bc_type_[face] = UNKNOWN;
  }


  return face_bc_[face];
}

BoundaryCond *GridInfo::copy_bc(BoundaryCond *bc, BoundaryCondition bc_type)
{
  BoundaryCond *ret = 0;

  switch (bc_type) {
  case SUBDOMAIN:
    ret = new SubdomainBc(dynamic_cast<const SubdomainBc&>(*bc));
    break;

  case EWALL:
    ret = new Ewall(dynamic_cast<const Ewall&>(*bc));
    break;

  case HWALL:
    ret = new Hwall(dynamic_cast<const Hwall&>(*bc));
    break;

  case PML:
    ret = new Pml(dynamic_cast<const Pml&>(*bc));
    break;

  case UNKNOWN:
    ret = new UnknownBc(dynamic_cast<const UnknownBc&>(*bc));
    break;
  }

  return ret;
}

unsigned int GridInfo::get_face_thickness(Face face)
{
  unsigned int ret = 0;

  if (face_bc_[face]) {
    ret = face_bc_[face]->get_thickness();
  }

  return ret;
}

void GridInfo::apply_boundaries(Grid &grid, FieldType type)
{
  for (unsigned int i = 0; i < 6; i++)
  {
    if (face_bc_[i]) {
      face_bc_[i]->apply(static_cast<Face>(i), grid, type);
    }
  }
}

