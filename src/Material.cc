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

#include "Material.hh"

Material::Material()
  : type_(NON_PERMEABLE), epsilon_(1), sigma_(0), 
    mu_(1), sigma_star_(0), name_("Free space"), vc_(0), fp_(0)
    // Look at me, I'm free space!
{
  
}

Material::~Material()
{
  
}

MaterialType Material::get_type() const
{
  return type_;
}

mat_prop_t Material::get_epsilon() const
{
  return epsilon_;
}

mat_prop_t Material::get_sigma() const
{
  return sigma_;
}

mat_prop_t Material::get_mu() const
{
  return mu_;
}

mat_prop_t Material::get_sigma_star() const
{
  return sigma_star_;
}

string Material::get_name() const
{
  return name_;
}

void Material::set_type(MaterialType type)
{
  type_ = type;
}

void Material::set_epsilon(mat_prop_t eps)
{
  epsilon_ = eps;
}

void Material::set_sigma(mat_prop_t sigma)
{
  sigma_ = sigma;
}

void Material::set_mu(mat_prop_t mu) 
{
  mu_ = mu;
}

void Material::set_sigma_star(mat_prop_t sigstar)
{
  sigma_star_ = sigstar;
}

void Material::set_name(string name)
{
  name_ = name;
}

void Material::set_collision_freq(mat_prop_t vc)
{
  vc_ = vc;
}

void Material::set_plasma_freq(mat_prop_t fp)
{
  fp_ = fp;
}

mat_prop_t Material::get_collision_freq() const
{
  return vc_;
}

mat_prop_t Material::get_plasma_freq() const
{
  return fp_;
}
