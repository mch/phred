
#include "SubdomainBc.hh"
#include <iostream>

using namespace std;

void SubdomainBc::apply(Face face, Grid &grid)
{
  cout << "Supposed to apply subdomain boundary condition; send interface to "
       << neighbour_ << endl;
}

void SubdomainBc::set_neighbour(int neighbour)
{
  neighbour_ = neighbour;
}

int SubdomainBc::get_neighbour()
{
  return neighbour_;
}
