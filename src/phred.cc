/* 
   phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004 Matt Hughes

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


/** \mainpage Phred Documentation
 *
 * \section Introduction
 *
 * Phred is a parallel finite difference time domain electromagnetics
 * simulator. It is written in C and C++, and uses MPI (Message
 * passing interface) for parallelism. It runs best on a cluster of
 * homogenous nodes, where the number of nodes is an even number. 
 *
 * \section addtocode Adding to the code
 * 
 * Phred is designed to be an extensible and flexible research
 * code. As such, it is hopefully easy for other programmers to add
 * boundary conditions, excitations, output formats, and other cool
 * things. This section summarizes how to add such items to Phred.
 *
 * \subsection addbc Adding a boundary condition
 *
 * \subsection adde Adding an excitation
 *
 * Phred has built in support for excitations which vary only with
 * time and for excitations which vary with time and
 * space. Excitations, by default, are defined over some rectangular
 * region of space. It's possible to change this by overriding the
 * Excite::excite() method of the Excitation class. 
 *
 * The classes TimeExcitation and SpaceExcitation take care of looping
 * over the default block of space by implementing the excite()
 * method. Sources which need to be applied to some other shape of
 * space can override Excitation::excite() to do so. 
 *
 * Any new excitation which is content to be applied to a block of the
 * grid needs to only override the source_function() method of
 * TimeExcitation or SpaceExcitation. 
 *
 * The source_function() method is public, so that other things, such
 * as the total/scattered thing can use the same excitation objects. 
 *
 *
 * \subsection addoutput Adding a output method
 *
 * Phred supports MATLAB friendly ASCII and binary, netCDF, and HDF
 * output formats. To add support for you favorite file format, make a
 * class which inherits from ... and implement the methods of the
 * interface. 
 */

/*#include <stdio.h>
#include <sys/types.h>
#include "system.h"*/

/* MPI (rocks your socks right off) */
#include <mpi.h>

#ifdef HAVE_LIBPOPT
/* popt plays way nicer with MPI than getopt. Trust me. */
#include <popt.h>
#endif

/* Let's use C++ for things that aren't speed critical, because life
   is just so much easier that way. And safer. Practice safe hex. */
#include <iostream>
#include <vector>
#include <string>

using namespace std; // Too lazy to type namespaces all the time. 

// In order to actually do stuff...
#include "MaterialLib.hh"
#include "Grid.hh"

#define EXIT_FAILURE 1

static void usage (int status);

#ifdef HAVE_LIBPOPT
static struct poptOption options[] = 
  {
    {"help", 'h', POPT_ARG_NONE, 0, 'h'},
    {"version", 'V', POPT_ARG_NONE, 0, 'V'},
    {"verbose", 'v', POPT_ARG_NONE, 0, 'v'},
    {"file", 'f', POPT_ARG_STRING, 0, 'f'}
  };
#endif

static int decode_switches (int argc, char **argv);

// Ugly globals
string inputfile;
const char *program_name;

int
main (int argc, char **argv)
{
  int i, rank, size, len;
  string prog_name;
  char *temp;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0)
  {
    prog_name = argv[0];
    len = prog_name.size();
    program_name = prog_name.c_str();
    
    // MPI implementations are not required to distribute command line
    // args, although MPICH does.
    i = decode_switches (argc, argv);
  } 

  // Our first try at MPI
  MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    temp = new char(len + 1);
  } else {
    temp = const_cast<char *>(prog_name.c_str());
  }

  MPI_Bcast(temp, len + 1, MPI_CHAR, 0, MPI_COMM_WORLD);

  if (rank != 0)
  {
    prog_name = temp;
    delete[] temp;
  }

  // Please sir, what file should I load? 
  if (rank == 0)
  {
    len = inputfile.size();
    temp = const_cast<char *>(inputfile.c_str());
  }

  MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    temp = new char(len + 1);
  }

  MPI_Bcast(temp, len + 1, MPI_CHAR, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    inputfile = temp;
    delete[] temp;
  }

  cout << "My rank is " << rank << ", and my program name is " 
       << prog_name << " which is " << len << " chars long." << endl;
  cout << "I'm going to load data from the file '" << inputfile << "'."
       << endl;

  // Parse the input script (each process will just load it's own file
  // for now. ) 

  // Subdomain the grid among the available processors and have each
  // processor set up its grid.



  // THIS STUFF HERE IS TEMPORARY
  Grid grid; 
  GridInfo info;
  
  info.global_dimx_ = info.dimx_ = 100;
  info.global_dimy_ = info.dimy_ = 100;
  info.global_dimz_ = info.dimz_ = 100;
  info.deltax_ = info.deltay_ = info.deltaz_ = 18.75e-9;
  info.deltat_ = 36e-18;
  info.start_x_ = info.start_y_ = info.start_z_ = 0;
  
  grid.setup_grid(info);
  grid.alloc_grid();

  MaterialLib mats; 
  Material mat; // defaults to free space
  mats.add_material(mat);
  
  // The library stores copies. 
  mat.set_epsilon(2.2);
  mat.set_name("Substrate");
  mats.add_material(mat);

  grid.load_materials(mats);

  // Global coordinates. 
  grid.define_box(0, 100, 0, 100, 0, 100, 1);
  grid.define_box(40, 60, 40, 60, 40, 60, 2);
  grid.set_boundary(FRONT, EWALL);
  grid.set_boundary(BACK, EWALL);
  grid.set_boundary(LEFT, EWALL);
  grid.set_boundary(RIGHT, EWALL);
  grid.set_boundary(BOTTOM, EWALL);
  grid.set_boundary(TOP, EWALL);

  // Main loop
  unsigned int num_time_steps = 100;
  unsigned int ts = 0;
  
  for (ts = 0; ts < num_time_steps; ts++) {
    // Excitations

    // E Field update

    // H Field update

    // Boundary condition application

    // Subdomain interface plane sharing

    // DFT outputs
    
    // Data block outputs
    
  }

  // Thank you and goodnight
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  exit (0);
}

/* Set all the option flags according to the switches specified.
   Return the index of the first non-option argument.  */

static int
decode_switches (int argc, char **argv)
{
  int c;
  char *arg = 0;

#ifdef HAVE_LIBPOPT
  poptContext ctx = poptGetContext(0, argc, 
                                   const_cast<const char **>(argv), 
                                   options, 0);

  while ((c = poptGetNextOpt (ctx)) > 0 || c == POPT_ERROR_BADOPT)
  {
    if (c == POPT_ERROR_BADOPT)
      continue;

    switch (c)
    {
    case 'V':
      cout << "phred " << VERSION << endl;
      exit (0);
      break;

    case 'h':
      usage (0);
      break;

    case 'f':
      arg = const_cast<char *>(poptGetOptArg(ctx));

      if (arg)
        inputfile = arg;
      else 
      {
        cout << "No filename given for the -f switch." << endl;
        exit(0);
      }
      break;

    default:
      cout << "WARNING: got unknown option number: " << c << endl;
    }
  }

  poptFreeContext(ctx);

#endif

  
  return optind;
}


static void
usage (int status)
{
  printf ("%s - \
Phred is a parallel finite difference time domain electromagnetics simulator.\n", program_name);

#ifdef HAVE_LIBPOPT
  printf ("Usage: %s [OPTION]... [FILE]...\n", program_name);
  printf ("\
Options:\n\
  -f, --filename             filename to read problem description from\n\
  -v, --verbose              print more information\n\
  -h, --help                 display this help and exit\n\
  -V, --version              output version information and exit\n\
"));
#else
  printf ("Usage: %s FILENAME\n", program_name);
  printf("\
FILENAME is the file containing the description of the proble to simulate.\n\
");
#endif

  exit (status);
}
