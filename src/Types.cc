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

/** \file Types.cc
 * Utility functions for operating on types defined in
 * Types.hh.in. Things like outputting to a stream, etc.
 */
 
#include "Types.hh"
#include <ostream>

using namespace std;

ostream& operator<< (ostream& os, const grid_point &p)
{
  return os << "(" << p.x << ", " << p.y << ", " << p.z << ")";
}

ostream& operator<< (ostream& os, const point &p)
{
  return os << "(" << p.x << ", " << p.y << ", " << p.z << ")";
}

const char *field_component_string(const FieldComponent &fc)
{
  switch(fc)
  {
  case FC_EX:
    return "Ex";
    break;

  case FC_EY:
    return "Ey";
    break;

  case FC_EZ:
    return "Ez";
    break;

  case FC_HX:
    return "Hx";
    break;

  case FC_HY:
    return "Hy";
    break;

  case FC_HZ:
    return "Hz";
    break;

  case FC_E:
    return "E intensity";
    break;

  case FC_H:
    return "H intensity";
    break;
  }

  return "invalid FieldComponent value";
}

const char *face_string(const Face &f)
{
  switch (f)
  {
  case FRONT:
    return "front";
    break;

  case BACK:
    return "back";
    break;

  case LEFT:
    return "left";
    break;

  case RIGHT:
    return "right";
    break;

  case BOTTOM:
    return "bottom";
    break;

  case TOP:
    return "top";
    break;
  }

  return "invalid Face value";
}
