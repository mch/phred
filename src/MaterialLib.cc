/* 
   phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004 Matt Hughes <mhughe@uvic.ca>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#include "MaterialLib.hh"
#include "Exceptions.hh"

MaterialLib::MaterialLib()
{
  Material temp;
  temp.set_epsilon(0);
  temp.set_sigma(INFINITY);
  temp.set_mu(0);

  add_material("PEC", temp);

  temp.set_epsilon(1);
  temp.set_sigma(0);
  temp.set_mu(1);
  
  add_material("FreeSpace", temp);
}

MaterialLib::~MaterialLib()
{
  
}

void MaterialLib::add_material(const char *name, Material &mat)
{
  mat.set_name(name);
  materials_[name] = mat;
}

int MaterialLib::num_materials() const
{
  return materials_.size();
}

const Material &MaterialLib::get_material(const char *name) const
{
  map<string, Material>::const_iterator iter = materials_.find(name);

  if (iter == materials_.end())
    throw UnknownMaterialException("Unknown material requested.");

  return (*iter).second;
}
