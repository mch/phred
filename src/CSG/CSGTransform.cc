/* 
   Phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004 Matt Hughes <mhughe@uvic.ca>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#include "CSGTransform.hh"
#include <math.h>

CSGTransform::CSGTransform(shared_ptr<CGSObject> child)
  : child_(child), tx_(0), ty_(0), tz_(0), angle_(0), 
    sx_(1), sy_(1), sz_(1)
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      A[i][j] = 0;
      Ainv[i][j] = 0;
    }
}

CSGTransform::~CSGTransform()
{}

bool CSGTransform::is_point_inside(float x, float y, float z) const
{
  // Transform the point to the child object's coordinates. 
  float x0 = x - tx_ + rotation_point_.x;
  float y0 = y - ty_ + rotation_point_.y;
  float z0 = z - tz_ + rotation_point_.z;

  float x1 = x0 * Ainv[0][0] + y0 * Ainv[0][1] + z0 * Ainv[0][2]
    - rotation_point_.x;
  float y1 = x0 * Ainv[1][0] + y0 * Ainv[1][1] + z0 * Ainv[1][2]
    - rotation_point_.y;
  float z1 = x0 * Ainv[2][0] + y0 * Ainv[2][1] + z0 * Ainv[2][2]
    - rotation_point_.z;

  return (*child_).is_point_inside(x1, y1, z1);
}

shared_ptr<CSGObject> CSGTransform::copy() const;
{
  return shared_ptr<CSGObject> (new CSGTransform(*this));
}

void CSGTransform::set_rotation(point p, point v, float angle)
{
  float norm = 1/sqrt(v.x * v.x + v.y * v.y + v.z * v.z);

  vector_.x = v.x * norm;
  vector_.y = v.y * norm;
  vector_.z = v.z * norm;

  rotation_point_ = p;
  angle_ = angle;

  // Calculate the transformation matrix. 
  calc_rotation_matrix(rotation_point_, vector_, angle_, &A);
  calc_rotation_matrix(rotation_point_, vector_, (-1) * angle_, &Ainv);
}

void CSGTransform::calc_rotation_matrix(const point &v, 
                                        const float &angle, 
                                        float **A) const
{
  float half_angle = angle / 2;
  float sha = sin(half_angle);
  float x = v.x * sha;
  float y = v.y * sha;
  float z = v.z * sha;
  float w = cos(half_angle);

  A[0][0] = 1 - 2 * (y*y + z*z);
  A[1][0] = 2 * (x*y + w*z);
  A[2][0] = 2 * (x*z - w*y);
  
  A[0][1] = 2 * (x*y - w*z);
  A[1][1] = 1 - 2*(x*x + z*z);
  A[2][1] = 2 * (y*z + w*x);

  A[0][2] = 2 * (x*z + w*y);
  A[1][2] = 2 * (y*z - w*x);
  A[2][2] = 1 - 2 * (x*x + y*y);
}

void set_scaling(float sx, float sy, float sz)
{
  sx_ = sx;
  sy_ = sy;
  sz_ = sz;
}

void set_translation(float tx, float ty, float tz)
{
  tx_ = tx;
  ty_ = ty;
  tz_ = tz;
}


