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

#ifndef CSG_OPERATOR_H
#define CSG_OPERATOR_H

#include "CSGObject.hh"
#include <boost/shared_ptr.hpp>

using namespace boost;

/**
 * This is an abstract base class for a constructive solid geometry
 * operator. Operators have two children, so this class has data
 * related to binary trees. 
 */ 
class CSGOperator : public CSGObject {
public:

  CSGOperator(shared_ptr<CSGObject> left, 
              shared_ptr<CSGObject> right)
    : left_(left), right_(right)
  {}

  virtual ~CSGOperator()
  {}

protected:
  shared_ptr<CSGObject> left_;
  shared_ptr<CSGObject> right_;
};

#endif // CSG_OBJECT_H
