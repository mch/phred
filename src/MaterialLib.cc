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

#include <math.h>
#include <fstream>

using namespace std;

MaterialLib::MaterialLib()
{
  Material temp, temp2;
  temp.set_epsilon(0);
  temp.set_sigma(INFINITY);
  temp.set_mu(0);
  add_material("PEC", temp);

  temp2.set_epsilon(1);
  temp2.set_sigma(0);
  temp2.set_mu(1);
  add_material("freespace", temp2);
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

void MaterialLib::load_material_file(const char *filename)
{
  
}

void MaterialLib::save_material_file(const char *filename)
{
  ofstream file;

  file.open(filename, ofstream::trunc | ofstream::out);
  
  if(file.is_open())
  {
    file << "Phred " << PACKAGE_VERSION << " Material Library\n\n";
    
    map<string, Material>::const_iterator iter = materials_.begin();
    map<string, Material>::const_iterator iter_e = materials_.end();

    for (; iter != iter_e; ++iter)
    {
      file << (*iter).second.get_name() << " " << (*iter).second.get_epsilon()
           << " " << (*iter).second.get_mu()
           << " " << (*iter).second.get_sigma()
           << " " << (*iter).second.get_sigma_star() << "\n";

      const map<string, mat_prop_t> &props = (*iter).second.get_properties();

      map<string, mat_prop_t>::const_iterator piter = props.begin();
      map<string, mat_prop_t>::const_iterator piter_e = props.end();

      for (; piter != piter_e; ++piter)
      {
        file << (*piter).first << " " << (*piter).second << "\n";
      }
      file << "\n";
    }

    file.close();
  } else {
    throw DataWriterException("Unable to open a file to write material database to.");
  }
}
