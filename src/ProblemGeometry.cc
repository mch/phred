/* 
   Phred - Phred is a parallel finite difference time domain
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

#include "ProblemGeometry.hh"
#include "Grid.hh"

ProblemGeometry::ProblemGeometry()
  : unit_(1), grid_box_(string("FreeSpace"), 
                        shared_ptr<CSGObject>(new CSGBox()))
{}

ProblemGeometry::~ProblemGeometry()
{}

void ProblemGeometry::set_grid_material(const char *material)
{
  grid_box_.material_ = material;
}

unsigned int 
ProblemGeometry::get_material_id(float x, float y, float z) const
{
  vector <GeomObject>::const_iterator iter = objects_.begin();
  vector <GeomObject>::const_iterator iter_e = objects_.end();
  
  unsigned int ret = grid_box_.material_id_;

  for(; iter != iter_e; ++iter)
  {
    if ((*((*iter).obj_)).is_point_inside(x, y, z) == INSIDE)
    {
      ret = (*iter).material_id_;
      break;
    }
  }

  return ret;
}

void ProblemGeometry::add_object(string material, shared_ptr<CSGObject> obj)
{
  objects_.push_back(GeomObject(material, obj));
}

void ProblemGeometry::set_grid_size(float x_size, float y_size, float z_size)
{
  CSGBox *box = dynamic_cast<CSGBox *>(grid_box_.obj_.get());

  if (box)
    box->set_size(x_size, y_size, z_size);
}

point ProblemGeometry::get_grid_size()
{
  CSGBox *box = dynamic_cast<CSGBox *>(grid_box_.obj_.get());
  point ret;

  if (box)
    box->get_size();

  return ret;
}

void ProblemGeometry::set_grid_centre(float x, float y, float z)
{
  CSGBox *box = dynamic_cast<CSGBox *>(grid_box_.obj_.get());

  if (box)
    box->set_centre(x, y, z);
}

point ProblemGeometry::get_grid_centre()
{
  CSGBox *box = dynamic_cast<CSGBox *>(grid_box_.obj_.get());
  point ret;

  if (box)
    ret = box->get_centre();

  return ret;
}

void ProblemGeometry::init(const Grid &grid)
{
  vector <GeomObject>::iterator iter = objects_.begin();
  vector <GeomObject>::iterator iter_e = objects_.end();
  const MaterialLib &material = grid.get_material_lib();

  try {
    const Material &mat = material.get_material(grid_box_.material_.c_str());
    grid_box_.material_id_ = mat.get_id();
  } catch (const UnknownMaterialException &e) {
    cout << "WARNING! The grid is using the material '"
         << (*iter).material_ << "' which does not exist!" << endl;
    throw e;
  }

  for(; iter != iter_e; ++iter)
  {
    try {
      const Material &mat = material.get_material((*iter).material_.c_str());

      (*iter).material_id_ = mat.get_id();

    } catch (const UnknownMaterialException &e) {
      cout << "WARNING! A solid object is using the material '"
           << (*iter).material_ << "' which does not exist!" << endl;
      throw e;
    }
  }  
}

void ProblemGeometry::deinit()
{

}
