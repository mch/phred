#include "MaterialLib.hh"

MaterialLib::MaterialLib()
{

}

MaterialLib::~MaterialLib()
{
  
}

void MaterialLib::add_material(Material &mat)
{
  materials_.push_back(mat);
}

int MaterialLib::num_materials()
{
  materials_.size();
}

vector<Material>::iterator MaterialLib::get_material_iter_begin()
{
  return materials_.begin();
}

vector<Material>::iterator MaterialLib::get_material_iter_end()
{
  return materials_.end();
}
