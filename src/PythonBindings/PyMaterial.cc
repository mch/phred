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

#include <boost/python.hpp>

#include "../MaterialLib.hh"

using namespace boost::python;

void export_materials()
{
  class_<Material>("Material", "Material class")
    .def("get_type", &Material::type)
    .add_property("epsilon", &Material::get_epsilon, &Material::set_epsilon)
    .add_property("sigma", &Material::get_sigma, &Material::set_sigma)
    .add_property("sigma_star", &Material::get_sigma_star, 
                  &Material::set_sigma_star)
    .add_property("mu", &Material::get_mu, &Material::set_mu)
    .add_property("name", &Material::get_name, &Material::set_name)
    .add_property("collision_freq", &Material::get_collision_freq, 
                  &Material::set_collision_freq)
    .add_property("plasma_freq", &Material::get_plasma_freq, 
                  &Material::set_plasma_freq)
    ;

  class_<MaterialLib>("MaterialLib", "Materials Library")
    .def("add_material", &MaterialLib::add_material)
    .def("num_materials", &MaterialLib::num_materials)
    .def("save_file", &MaterialLib::save_material_file)
    .def("load_file", &MaterialLib::load_material_file)
    ;
}
