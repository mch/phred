#include "Material.hh"

Material::Material()
  : type_(NON_PERMEABLE), epsilon_(1), sigma_(0), 
    mu_(1), sigma_star_(0), name_("Free space")
    // Look at me, I'm free space!
{
  
}

Material::~Material()
{
  
}

MaterialType Material::get_type()
{
  return type_;
}

mat_prop_t Material::get_epsilon()
{
  return epsilon_;
}

mat_prop_t Material::get_sigma()
{
  return sigma_;
}

mat_prop_t Material::get_mu()
{
  return mu_;
}

mat_prop_t Material::get_sigma_star()
{
  return sigma_star_;
}

string Material::get_name()
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
