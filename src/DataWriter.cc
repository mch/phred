#include "DataWriter.hh"

Variable_4d& DatWriter::add_4d_variable(string name, region_t r)
{
  if (variables_.size() >= max_vars_)
    throw std::exception("This DataWriter cannot support more variables.");
  
  Variable_4d *var = new Variable_4d(*this, name, r);
  
  // BUG: Check for duplicates!!!!
  variables_[name] = r;

  return *var;
}


void Variable_4d::save_point(unsigned int x, unsigned int y,
                             unsigned int z, field_t val)
{
  if (dw_.get_rank() > 0) {

  } else {

  }
}
