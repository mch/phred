#include "Gaussm.hh"

Gaussm::Gaussm()
  : alpha_(1), deltaf_(50e6), f0_(1e9)
{}

Gaussm::~Gaussm()
{}

field_t Gaussm::source_function(Grid &grid, unsigned int time_step)
{
  field_t t = time_step * grid.get_deltat();

  return alpha_ * exp(-pow((t - 4. / (PI * deltaf_)) * deltaf_ * PI, 2))
    * sin(2. * PI * f0_ * (t - 4. / (PI * deltaf_)));
}

void Gaussm::set_parameters(field_t alpha, field_t deltaf, field_t f0)
{
  alpha_ = alpha;
  deltaf_ = deltaf;
  f0_ = f0;
}

field_t Gaussm::get_alpha()
{
  return alpha_;
}

field_t Gaussm::get_deltaf()
{
  return deltaf_;
}

field_t Gaussm::get_f0()
{
  return f0_;
}


