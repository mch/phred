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

#include "PyInterpreter.hh"
#include "../Exceptions.hh"

#include "PyTypes.hh"
#include "PyMaterial.hh"
#include "PyResults.hh"
#include "PyGeometry.hh"
#include "PyExcitations.hh"
#include "PyFDTD.hh"
#include "PyDatawriters.hh"
#include "PyBoundaries.hh"
//#include "PyGrid.hh"

static struct _inittab modules_[] = 
  {
    {"Results", &initResults},
    {"Excitations", &initExcitations},
    {"FDTD", &initFDTD},
    {"Boundaries", &initBoundaries},
    {"Materials", &initMaterials},
    {"Types", &initTypes},
    {"DataWriters", &initDataWriters},
    {"Geometry", &initGeometry},
    //{"Grid", &initGrid},
    {0, 0}
  };
  
#include <stdio.h>

#ifdef HAVE_LIBREADLINE
#include <unistd.h>
#include <readline/readline.h>
#include <readline/history.h>
#endif

#include <string>

using namespace std;

PyInterpreter::PyInterpreter(int rank, int size)
  : rank_(rank), size_(size)
{
  int ret = PyImport_ExtendInittab(modules_);
  if (ret == -1)
    throw PyInterpException("Unable to add Phred modules to Python interpreter.");
  
  Py_Initialize();

#ifdef HAVE_LIBREADLINE
  rl_bind_key ('\t', rl_insert);
#endif
}

PyInterpreter::~PyInterpreter()
{
  Py_Finalize();
}

void PyInterpreter::run_script(const char *filename)
{
  // Have to use stdio because I need the FILE ptr for Python.
  FILE *fp;
  fp = fopen(filename, "r");

  if (fp == 0)
    throw PyInterpException("Unable to open Python script file.");

  add_modules();

  PyRun_SimpleFile(fp, filename);
  fclose(fp);
}

void PyInterpreter::run()
{
  if (rank_ == 0)
    master();
  else
    slave();
}

void PyInterpreter::slave()
{
  // Recieving a message of zero length means that it is time to
  // return. -1 means that the root rank is not bound to a terminal
  // and that we should throw one too. 
  int size = 1;
  char *buffer = 0;

  handle<> main_module(borrowed( PyImport_AddModule("__main__") ));
  handle<> main_namespace(borrowed( PyModule_GetDict(main_module.get()) ));
  
  add_modules();

  while (size > 0) 
  {
    MPI_Bcast(static_cast<void *>(&size), 1, MPI_INT, 
              0, MPI_COMM_WORLD);

    if (size == -1)
      throw PyInterpException("Standard input and standard output must be bound to a terminal to use interactive mode.");

    // A wasty memory management policy. But it's simple and should
    // always work, barring lack of memory.
    buffer = new char[size + 1];

    if (!buffer)
      throw MemoryException();

    MPI_Bcast(static_cast<void *>(buffer), size, MPI_CHAR, 
              0, MPI_COMM_WORLD);
    buffer[size] = 0;

    try {
      handle<> res( PyRun_String(buffer, Py_file_input, 
                                 main_namespace.get(),
                                 main_namespace.get()));
    }
    catch (error_already_set)
    {
      PyErr_Print();
    }

    delete[] buffer;
  }
}

void PyInterpreter::master()
{
  // If stdin isn't a tty, don't even try. 
  if (!isatty(STDIN_FILENO)) {
    // Tell slaves to abort as well. 
    int size = -1;
    MPI_Bcast(static_cast<void *>(&size), 1, MPI_INT, 0, MPI_COMM_WORLD);    

    throw PyInterpException("Standard input and standard output must be bound to a terminal to use interactive mode.");
  }

  handle<> main_module(borrowed( PyImport_AddModule("__main__") ));
  handle<> main_namespace(borrowed( PyModule_GetDict(main_module.get()) ));
  
  add_modules();

  //if (size_ > 1) {
#ifdef HAVE_LIBREADLINE
    cout << "Phred interactive Python interpreter running. Type ctrl-d to quit." << endl;
    char *ln = 0; 
    const char *prompt, *p1 = ">>> ", *p2 = "... ";
    prompt = p1;
    bool multiline = false;
    int slen = 0;
    string buffer;
    while ((ln = rl())) {
      buffer += ln;
      buffer += '\n';

      // Ugg, need to check if the last non white space char is a ':', 
      // if so, go multiline, also need to ignore comments. 
      if (strlen(ln) > 0 && ln[strlen(ln) - 1] == ':') 
      {
        prompt = p2;
        multiline = true;
      }

      if ((multiline && strlen(ln) == 0) || (!multiline && buffer.size() > 0))
      {
        int size = buffer.length();
        MPI_Bcast(static_cast<void *>(&size), 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(static_cast<void *>(const_cast<char *>(buffer.c_str())), buffer.length(), 
                  MPI_CHAR, 0, MPI_COMM_WORLD);

        try {
          handle<> res( PyRun_String(buffer.c_str(), Py_single_input, 
                                     main_namespace.get(),
                                     main_namespace.get()));
        }
        catch (error_already_set)
        {
          PyErr_Print();
        }
        multiline = false;
        prompt = p1;
        buffer.clear();
      }
    }
    cout << endl;
#else
    throw PyInterpException("Only one process can be used in interactive mode.");
#endif
    //}  
  // else
//   {
//     char **argv;
//     argv = new char*[2];
//     argv[0] = "phred";
//     argv[1] = "\0";

//     Py_Main(1, argv);  

//     delete[] argv;
//   }
}

void PyInterpreter::add_modules()
{
  handle<> main_module(borrowed( PyImport_AddModule("__main__") ));
  handle<> main_namespace(borrowed( PyModule_GetDict(main_module.get()) ));
  
  handle<> rname ( PyString_FromString("Results") );
  handle<> res( PyImport_Import(rname.get()) );
  PyDict_SetItemString(main_namespace.get(), "Results", res.get());

  handle<> ename ( PyString_FromString("Excitations") );
  handle<> ex( PyImport_Import(ename.get()) );
  PyDict_SetItemString(main_namespace.get(), "Excitations", ex.get());

  handle<> fname ( PyString_FromString("FDTD") );
  handle<> fdtd( PyImport_Import(fname.get()) );
  PyDict_SetItemString(main_namespace.get(), "FDTD", fdtd.get());

  handle<> mname ( PyString_FromString("Materials") );
  handle<> mats( PyImport_Import(mname.get()) );
  PyDict_SetItemString(main_namespace.get(), "Materials", mats.get());

  handle<> tname ( PyString_FromString("Types") );
  handle<> types( PyImport_Import(tname.get()) );
  PyDict_SetItemString(main_namespace.get(), "Types", types.get());

  handle<> bname ( PyString_FromString("Boundaries") );
  handle<> bnds( PyImport_Import(bname.get()) );
  PyDict_SetItemString(main_namespace.get(), "Boundaries", bnds.get());

  handle<> gname ( PyString_FromString("Geometry") );
  handle<> geom( PyImport_Import(gname.get()) );
  PyDict_SetItemString(main_namespace.get(), "Geometry", geom.get());

  handle<> dwname ( PyString_FromString("DataWriters") );
  handle<> dws( PyImport_Import(dwname.get()) );
  PyDict_SetItemString(main_namespace.get(), "DataWriters", dws.get());

  // MPI Data
  PyDict_SetItemString(main_namespace.get(), "MPI_RANK", 
                       PyInt_FromLong(rank_));
  PyDict_SetItemString(main_namespace.get(), "MPI_SIZE", 
                       PyInt_FromLong(size_));
}

char *PyInterpreter::rl()
{
#ifdef HAVE_LIBREADLINE
  char *line_read = readline (">> ");

  /* If the line has any text in it,
     save it on the history. */
  if (line_read && *line_read)
    add_history (line_read);

  return (line_read);

#else
  return 0;
#endif
}
