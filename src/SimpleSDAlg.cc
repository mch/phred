#include "SimpleSDAlg.hh"

SimpleSDAlg::SimpleSDAlg()
{}

SimpleSDAlg::~SimpleSDAlg()
{}

GridInfo SimpleSDAlg::decompose_domain(int rank, int size, 
                                   const GridInfo &info)
{
  bool divided = false;
  
  // n, m, and p are the number of divisions along the x, y, and z
  // axis' respectivly
  unsigned int sdx, sdy, sdz, n, m, p;
  sdx = info.global_dimx_;
  sdy = info.global_dimy_;
  sdz = info.global_dimz_;

  n = m = p = 1;

  for (int i = 0; i < size; i++)
  {
    if (sdx >= sdy && sdx >= sdz && (n+1)*m*p <= size) {
      n++;
      sdx = info.global_dimx_ / n;
      divided = true;
    }

    else if (sdy >= sdx && sdy >= sdz && n*(m+1)*p <= size) {
      m++;
      sdy = info.global_dimy_ / m;
      divided = true;
    }

    else if (sdz >= sdx && sdz >= sdy && n*m*(p+1) <= size) {
      p++;
      sdz = info.global_dimz_ / p;
      divided = true;
    }

    if (!divided) {
      if ((n+1)*m*p <= size) {
        n++;
        sdx = info.global_dimx_ / n;
      }
      else if (n*(m+1)*p <= size) {
        m++;
        sdy = info.global_dimy_ / m;
      }
      else if (n*m*(p+1) <= size) {
        p++;
        sdz = info.global_dimz_ / p;
      }
    }

    divided = false;
  }

  if (n*m*p != size) {
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

  // Calculate the starting point and length of the axis' including
  // the overlap required between subdomains. Assign boundary
  // conditions to the grid. 
  
}
