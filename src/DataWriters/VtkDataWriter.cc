/* 
   phred - Phred is a parallel finite difference time domain
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

#include "VtkDataWriter.hh"

#ifdef USE_VTK

VtkDataWriter::VtkDataWriter(int rank, int size)
  : DataWriter(rank, size)
{
  writer_.SetDataModeToBinary();
}

VtkDataWriter::~VtkDataWriter()
{

}

VtkDataWriter::init(const Grid &grid)
{

}

VtkDataWriter::deinit(const Grid &grid)
{

}

VtkDataWriter::add_variable()
{
  
}

unsigned int VtkDataWriter::write_data(unsigned int time_step, 
                                       Data &data, MPI_Datatype t, 
                                       void *ptr, unsigned int len)
{
  
}

#endif 