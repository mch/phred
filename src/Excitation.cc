#include "Excitation.hh"

Excitation::Excitation()
  : component_(EY)
{}

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
