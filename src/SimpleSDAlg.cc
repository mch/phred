#include "SimpleSDAlg.hh"
#include <exception>

SimpleSDAlg::SimpleSDAlg()
{}

SimpleSDAlg::~SimpleSDAlg()
{}

GridInfo SimpleSDAlg::decompose_domain(int rank, int size, 
                                       GridInfo &info)
{
  bool divided = false;

  if (size < 0) // that just wrong
    throw std::exception();

  unsigned int sz = static_cast<unsigned int>(size);
  
  // n, m, and p are the number of divisions along the x, y, and z
  // axis' respectivly
  unsigned int sdx, sdy, sdz, n, m, p;
  sdx = info.global_dimx_;
  sdy = info.global_dimy_;
  sdz = info.global_dimz_;

  n = m = p = 1;

  for (int i = 0; i < size; i++)
  {
    if (sdx >= sdy && sdx >= sdz && (n+1)*m*p <= sz) {
      n++;
      sdx = info.global_dimx_ / n;
      divided = true;
    }

    else if (sdy >= sdx && sdy >= sdz && n*(m+1)*p <= sz) {
      m++;
      sdy = info.global_dimy_ / m;
      divided = true;
    }

    else if (sdz >= sdx && sdz >= sdy && n*m*(p+1) <= sz) {
      p++;
      sdz = info.global_dimz_ / p;
      divided = true;
    }

    if (!divided) {
      if ((n+1)*m*p <= sz) {
        n++;
        sdx = info.global_dimx_ / n;
      }
      else if (n*(m+1)*p <= sz) {
        m++;
        sdy = info.global_dimy_ / m;
      }
      else if (n*m*(p+1) <= sz) {
        p++;
        sdz = info.global_dimz_ / p;
      }
    }

    divided = false;
  }

  if (n*m*p != sz) {
    cerr << "WARNING: simple domain decomposition did not result in one grid per\nprocessor!" << endl;
  }

  GridInfo result;
  result.global_dimx_ = info.global_dimx_;
  result.global_dimy_ = info.global_dimy_;
  result.global_dimz_ = info.global_dimz_;
  
  result.deltax_ = info.deltax_;
  result.deltay_ = info.deltay_;
  result.deltaz_ = info.deltaz_;
  result.deltat_ = info.deltat_;


  // Find the domain this grid will be in. 
  unsigned int x, y, z; 
  x = rank % n;
  y = ((rank - x)/n) % m;
  z = ((rank - x)/n - y) / m;

  // Assign sizes and starting points including overlap
  result.dimx_ = static_cast<unsigned int>
    (floor(static_cast<double>(result.global_dimx_) / n));
  result.dimy_ = static_cast<unsigned int>
    (floor(static_cast<double>(result.global_dimy_) / m));
  result.dimz_ = static_cast<unsigned int>
    (floor(static_cast<double>(result.global_dimz_) / p));

  if (x == n-1) { // in case the number of cells isn't evenly divisible
    result.dimx_ = result.global_dimx_ - (n - 1)
      * static_cast<unsigned int>
      (floor(static_cast<double>(result.global_dimx_) / n));
  } 
  
  if (y == m-1) {
    result.dimy_ = result.global_dimy_ - (m - 1)
      * static_cast<unsigned int>
      (floor(static_cast<double>(result.global_dimy_) / m));
  } 

  if (z == p-1) {
    result.dimz_ = result.global_dimz_ - (p - 1)
      * static_cast<unsigned int>
      (floor(static_cast<double>(result.global_dimz_) / p));
  } 

  result.start_x_ = x * result.dimx_;
  result.start_y_ = y * result.dimx_;
  result.start_z_ = z * result.dimx_;
  
  // Assign boundary conditions the ranks to talk to 
  SubdomainBc *sdbc = 0;

  if (x == 0) { // BACK
    result.set_boundary(BACK, info.get_bc_type(BACK));
  } else {
    sdbc = &(dynamic_cast<SubdomainBc &>(result.set_boundary(BACK, 
                                                             SUBDOMAIN)));
    sdbc->set_neighbour((z*m + y) * n + (x-1));
    sdbc->set_rank(rank);
    result.dimx_++;
    result.start_x_--;
  }

  if (x == n - 1) { // FRONT
    result.set_boundary(FRONT, info.get_bc_type(FRONT));
  } else {
    sdbc = &(dynamic_cast<SubdomainBc &>(result.set_boundary(FRONT, 
                                                             SUBDOMAIN)));
    sdbc->set_neighbour((z*m + y) * n + (x+1));
    sdbc->set_rank(rank);
    result.dimx_++;
  } 

  if (y == 0) { // LEFT 
    result.set_boundary(LEFT, info.get_bc_type(LEFT));
  } else {
    sdbc = &(dynamic_cast<SubdomainBc &>(result.set_boundary(LEFT, 
                                                             SUBDOMAIN)));
    sdbc->set_neighbour((z*m + (y-1)) * n + x);
    sdbc->set_rank(rank);
    result.dimy_++;
    result.start_y_--;
  } 

  if (y == m - 1) { // RIGHT
    result.set_boundary(RIGHT, info.get_bc_type(RIGHT));
  } else {
    sdbc = &(dynamic_cast<SubdomainBc &>(result.set_boundary(RIGHT, 
                                                             SUBDOMAIN)));
    sdbc->set_neighbour((z*m + (y+1)) * n + x);
    sdbc->set_rank(rank);
    result.dimy_++;
  }

  if (z == 0) { // BOTTOM
    result.set_boundary(BOTTOM, info.get_bc_type(BOTTOM));
  } else {
    sdbc = &(dynamic_cast<SubdomainBc &>(result.set_boundary(BOTTOM, 
                                                             SUBDOMAIN)));
    sdbc->set_neighbour(((z-1)*m + y) * n + x);
    sdbc->set_rank(rank);
    result.dimz_++;
    result.start_z_--;
  }

  if (z == p - 1) { // TOP
    result.set_boundary(TOP, info.get_bc_type(TOP));
  } else {
    sdbc = &(dynamic_cast<SubdomainBc &>(result.set_boundary(TOP, 
                                                             SUBDOMAIN)));
    sdbc->set_neighbour(((z+1)*m + y) * n + x);
    sdbc->set_rank(rank);
    result.dimz_++;
  }

  return result;
}
