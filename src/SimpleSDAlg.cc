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
    result.dimx_ = result.global_dimx_ - n 
      * static_cast<unsigned int>
      (floor(static_cast<double>(result.global_dimx_) / n));
  } 
  
  if (y == m-1) {
    result.dimy_ = result.global_dimy_ - m 
      * static_cast<unsigned int>
      (floor(static_cast<double>(result.global_dimy_) / m));
  } 

  if (z == p-1) {
    result.dimz_ = result.global_dimz_ - p 
      * static_cast<unsigned int>
      (floor(static_cast<double>(result.global_dimz_) / p));
  } 

  result.start_x_ = x * result.dimx_;
  result.start_y_ = y * result.dimx_;
  result.start_z_ = z * result.dimx_;
  
  // Assign boundary conditions the ranks to talk to 
  if (x == 0) { // BACK
    result.face_bc_[BACK] = info.face_bc_[BACK];
  } else {
    result.face_bc_[BACK] = SUBDOMAIN;
    result.face_rank_[BACK] = (z*m + y) * n + (x-1);
    result.dimx_++;
    result.start_x_--;
  }

  if (x == n - 1) { // FRONT
    result.face_bc_[FRONT] = info.face_bc_[FRONT];    
  } else {
    result.face_bc_[FRONT] = SUBDOMAIN;
    result.face_rank_[FRONT] = (z*m + y) * n + (x+1);
    result.dimx_++;
  } 

  if (y == 0) { // LEFT 
    result.face_bc_[LEFT] = info.face_bc_[LEFT];
  } else {
    result.face_bc_[LEFT] = SUBDOMAIN;
    result.face_rank_[LEFT] = (z*m + (y-1)) * n + x;
    result.dimy_++;
    result.start_y_--;
  } 

  if (y == m - 1) { // RIGHT
    result.face_bc_[RIGHT] = info.face_bc_[RIGHT];
  } else {
    result.face_bc_[RIGHT] = SUBDOMAIN;
    result.face_rank_[RIGHT] = (z*m + (y+1)) * n + x;
    result.dimy_++;
  }

  if (z == 0) { // BOTTOM
    result.face_bc_[BOTTOM] = info.face_bc_[BOTTOM];
  } else {
    result.face_bc_[BOTTOM] = SUBDOMAIN;
    result.face_rank_[BOTTOM] = ((z-1)*m + y) * n + x;
    result.dimz_++;
    result.start_z_--;
  }

  if (z == p - 1) { // TOP
    result.face_bc_[TOP] = info.face_bc_[TOP];
  } else {
    result.face_bc_[TOP] = SUBDOMAIN;
    result.face_rank_[TOP] = ((z+1)*m + y) * n + x;
    result.dimz_++;
  }

  return result;
}
