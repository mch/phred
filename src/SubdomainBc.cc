
#include "SubdomainBc.hh"

void SubdomainBc::apply(Face face, Grid &grid)
{

}

void SubdomainBc::set_neighbour(int neighbour)
{
  neighbour_ = neighbour;
}

int SubdomainBc::get_neighbour()
{
  return neighbour_;
}
