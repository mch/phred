/* 
   Phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004-2005 Matt Hughes <mhughe@uvic.ca>

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

#include "Periodic.hh"

Periodic::Periodic(shared_ptr<PeriodicExcitation> pe)
  : pe_(pe), valid_(false)
{
  for (int i = 0; i < 6; i++)
  {
    faces_[i] = false;
  }
}

Periodic::~Periodic()
{

}

void Periodic::apply(Face face, Grid &grid, FieldType type)
{
  if (!valid_)
    throw BoundaryConditionException("Periodic: boundaries are not valid.");

  if (faces_[face])
  {
    if (type == E && (face == BACK || face == LEFT || face == BOTTOM))
      copy_e(face, grid);

    if (type == H && (face == TOP || face == RIGHT || face == FRONT))
      copy_h(face, grid);
  }
}

void Periodic::init(const Grid &grid, Face face)
{
  // This is where we find out what planes we've been applied to. 
  faces_[face] = true;

  bool fbv = false;
  bool lrv = false;
  bool tbv = false; 

  if (faces_[FRONT] && faces_[BACK]) 
    fbv = true; 

  if (faces_[LEFT] && faces_[RIGHT])
    lrv = true;

  if (faces_[TOP] && faces_[BOTTOM])
    tbv = true; 

  Face ex_face = pe_->get_face();

  if ((ex_face == FRONT || ex_face == BACK) && fbv)
    valid_ = false;

  else if ((ex_face == LEFT || ex_face == RIGHT) && lrv)
    valid_ = false;

  else if ((ex_face == TOP || ex_face == BOTTOM) && tbv)
    valid_ = false;

  else if (fbv || tbv || lrv)
    valid_ = true; 
}

void Periodic::deinit(const Grid &grid, Face face)
{

}

void Periodic::copy_e(Face face, Grid &grid)
{
  switch (face)
  {
  case BACK:
    for (int j = 0; j < grid.get_ldy_sd(); j++)
    {
      for (int k = 0; k < grid.get_ldz_sd(); k++)
      {
        grid.set_ey(0, j, k, grid.get_ey(grid.get_ldx_sd() - 1, j, k)); 
        grid.set_ez(0, j, k, grid.get_ez(grid.get_ldx_sd() - 1, j, k));
      }
    }

    

    break;

  case LEFT:
    for (int i = 0; i < grid.get_ldx_sd(); i++)
    {
      for (int k = 0; k < grid.get_ldz_sd(); k++)
      {
        grid.set_ex(i, 0, k, grid.get_ex(i, grid.get_ldy_sd() - 1, k));
        grid.set_ez(i, 0, k, grid.get_ez(i, grid.get_ldy_sd() - 1, k));
      }
    }
    break;

  case BOTTOM:
    for (int i = 0; i < grid.get_ldx_sd(); i++)
    {
      for (int j = 0; j < grid.get_ldy_sd(); j++)
      {
        grid.set_ex(i, j, 0, grid.get_ex(i, j, grid.get_ldz_sd() - 1));
        grid.set_ey(i, j, 0, grid.get_ey(i, j, grid.get_ldz_sd() - 1));
      }
    }
    break;
  }
  
}


void Periodic::copy_h(Face face, Grid &grid)
{
  switch (face)
  {
  case FRONT:
    for (int j = 0; j < grid.get_ldy_sd(); j++)
    {
      for (int k = 0; k < grid.get_ldz_sd(); k++)
      {
        grid.set_hy(grid.get_ldx_sd() - 1, j, k, 
                    grid.get_hy(0, j, k));

        grid.set_hz(grid.get_ldx_sd() - 1, j, k, 
                    grid.get_hz(0, j, k));
      }
    }
    break;

  case RIGHT:
    for (int i = 0; i < grid.get_ldx_sd(); i++)
    {
      for (int k = 0; k < grid.get_ldz_sd(); k++)
      {
        grid.set_hx(i, grid.get_ldy_sd() - 1, k, 
                    grid.get_hx(i, 0, k));
        grid.set_hz(i, grid.get_ldy_sd() - 1, k, 
                    grid.get_hz(i, 0, k));
      }
    }
    break;

  case TOP:
    for (int i = 0; i < grid.get_ldx_sd(); i++)
    {
      for (int j = 0; j < grid.get_ldy_sd(); j++)
      {
        grid.set_hx(i, j, grid.get_ldz_sd() - 1, 
                    grid.get_hx(i, j, 0));
        grid.set_hy(i, j, grid.get_ldz_sd() - 1, 
                    grid.get_hy(i, j, 0));
      }
    }
    break;

  }

}
