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

  // Result initially has the exact same properties as the input.
  GridInfo result = info; 

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

  if (x != 0) { // BACK
    sdbc = new SubdomainBc();

    sdbc->set_neighbour((z*m + y) * n + (x-1));
    sdbc->set_rank(rank);
    result.set_boundary(BACK, sdbc, true);

    result.dimx_++;
    result.start_x_--;
  }

  if (x != n - 1) { // FRONT
    sdbc = new SubdomainBc();

    sdbc->set_neighbour((z*m + y) * n + (x+1));
    sdbc->set_rank(rank);
    result.set_boundary(FRONT, sdbc, true);
    result.dimx_++;
  } 

  if (y != 0) { // LEFT 
    sdbc = new SubdomainBc();

    sdbc->set_neighbour((z*m + (y-1)) * n + x);
    sdbc->set_rank(rank);
    result.set_boundary(LEFT, sdbc, true);
    result.dimy_++;
    result.start_y_--;
  } 

  if (y != m - 1) { // RIGHT
    sdbc = new SubdomainBc();

    sdbc->set_neighbour((z*m + (y+1)) * n + x);
    sdbc->set_rank(rank);
    result.set_boundary(RIGHT, sdbc, true);
    result.dimy_++;
  }

  if (z != 0) { // BOTTOM
    sdbc = new SubdomainBc();

    sdbc->set_neighbour(((z-1)*m + y) * n + x);
    sdbc->set_rank(rank);
    result.set_boundary(BOTTOM, sdbc, true);
    result.dimz_++;
    result.start_z_--;
  }

  if (z != p - 1) { // TOP
    sdbc = new SubdomainBc();

    sdbc->set_neighbour(((z+1)*m + y) * n + x);
    sdbc->set_rank(rank);
    result.set_boundary(TOP, sdbc, true);
    result.dimz_++;
  }

  return result;
}
