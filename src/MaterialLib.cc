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
  return materials_.size();
}

vector<Material>::iterator MaterialLib::get_material_iter_begin()
{
  return materials_.begin();
}

vector<Material>::iterator MaterialLib::get_material_iter_end()
{
  return materials_.end();
}
