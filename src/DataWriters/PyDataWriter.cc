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

#include "PyDataWriter.hh"

PyDataWriter::PyDataWriter(int rank, int size)
  : DataWriter(rank, size)
{
#ifdef USE_PY_BINDINGS
  // Ensure that this class is only used if the script being run is a
  // Python script, and not a Jan script.
  if (!Py_IsInitialized())
    throw NoPythonException();
  
  

#else
    throw NoPythonException();
#endif
}

PyDataWriter::~PyDataWriter()
{

}

#ifdef USE_PY_BINDINGS
void PyDataWriter::init(const Grid &grid)
{}

void PyDataWriter::deinit(const Grid &grid)
{}

void PyDataWriter::add_variable(Result &result)
{
  const map<string, Variable *> vars = result.get_variables();
  map<string, Variable *>::const_iterator iter;
  map<string, Variable *>::const_iterator iter_e = vars.end();

  handle<> main_module( PyModule_New("__main__") );
  handle<> main_namespace(borrowed( PyModule_GetDict(main_module.get()) ));

  for (iter = vars.begin(); iter != iter_e; ++iter)
  {
    Variable *var = iter->second;
    const vector<Dimension> &dimensions = var->get_dimensions();
    
    if (dimensions.size() == 0)
      throw DataWriterException("Result must have at least one dimension.");
    
    numeric::array arr(make_tuple(1));
    vars_[var->get_name()] = arr;

    PyDict_SetItemString(main_namespace.get(), 
                         var->get_name().c_str(), arr.ptr());

  }  
}

unsigned int PyDataWriter::write_data(unsigned int time_step, 
                                      Variable &var, MPI_Datatype t, 
                                      void *ptr, unsigned int len)
{

}

#endif
