#include "Excitation.hh"

Excitation::Excitation()
  : x_start_(0), y_start_(0), z_start_(0), 
    x_end_(0), y_end_(0), z_end_(0), type_(E)
{
  polarization_[0] = 1;
  polarization_[1] = 0;
  polarization_[2] = 0;
}

Excitation::~Excitation()
{}

void Excitation::set_region(unsigned int x_start, unsigned int x_end, 
                            unsigned int y_start, unsigned int y_end, 
                            unsigned int z_start, unsigned int z_end)
{
  x_start_ = x_start;
  y_start_ = y_start;
  z_start_ = z_start;

  x_end_ = x_end;
  y_end_ = y_end;
  z_end_ = z_end;
}

void Excitation::set_polarization(field_t x, field_t y, field_t z)
{
  polarization_[0] = x;
  polarization_[1] = y;
  polarization_[2] = z;
}

void Excitation::set_type(FieldType t)
{
  type_ = t;
}
