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

#include "CSGOperator.hh"


CSGOperator::CSGOperator(shared_ptr<CSGObject> left, 
                         shared_ptr<CSGObject> right)
  : left_(left), right_(right)
{}

CSGOperator::~CSGOperator()
{}

CSGOperator::CSGOperator(const CSGOperator &rhs)
{
  *this = rhs;
}

const CSGOperator &CSGOperator::operator=(const CSGOperator &rhs)
{
  left_ = (*rhs.left_).copy();
  right_ = (*rhs.right_).copy();
  
  return *this;
}

shared_ptr<CSGObject> CSGOperator::copy() const
{
  CSGOperator *op = new CSGOperator((*left_).copy(), (*right_).copy());
  return shared_ptr<CSGOperator>(op);
}

std::ostream& CSGOperator::to_string(std::ostream &os) const
{
	return os << "CSGOperator...";
}

