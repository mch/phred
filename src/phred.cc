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
 * things. This section summarizes how to add such items to Phred in
 * C++. Phred embeds the Python language, and such items may also be
 * prototyped in Python by creating derived classes, although for
 * speed the implementation should be done in C/C++. 
 *
 * \subsection addbc Adding a boundary condition
 *
 * Boundary conditions can be added by subclassing the
 * BoundaryCondition object. For boundary conditions which only need
 * to consider the tangential and normal field components at a face,
 * adapter classes can be used to simplify the expression of the
 * boundary condition. See the electric (Ewall) and magnetic wall
 * boundary conditions (Hwall) for details. 
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

/* MPI (rocks your socks right off) */
#include <mpi.h>

/* Let's use C++ for things that aren't speed critical, because life
   is just so much easier that way. And safer. Practice safe hex. */
#include <iostream>
#include <vector>
#include <string>

using namespace std; // Too lazy to type namespaces all the time. 

// In order to actually do stuff...
#include "MaterialLib.hh"
#include "Grid.hh"
#include "SimpleSDAlg.hh"
#include "Gaussm.hh"
#include "config.h"
#include "Ewall.hh"
#include "Hwall.hh"
#include "PointResult.hh"
#include "AsciiDataWriter.hh"


#define EXIT_FAILURE 1

static void usage (int status);

#ifdef HAVE_LIBPOPT

/* popt plays way nicer with MPI than getopt. Trust me. */
#include <popt.h>

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
//   MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);

//   if (rank != 0) {
//     temp = new char(len + 1);
//   } else {
//     temp = const_cast<char *>(prog_name.c_str());
//   }

//   MPI_Bcast(temp, len + 1, MPI_CHAR, 0, MPI_COMM_WORLD);

//   if (rank != 0)
//   {
//     prog_name = temp;
//     delete[] temp;
//   }

//   // Please sir, what file should I load? 
//   if (rank == 0)
//   {
//     len = inputfile.size();
//     temp = const_cast<char *>(inputfile.c_str());
//   }

//   MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);

//   if (rank != 0) {
//     temp = new char(len + 1);
//   }

//   MPI_Bcast(temp, len + 1, MPI_CHAR, 0, MPI_COMM_WORLD);

//   if (rank != 0) {
//     inputfile = temp;
//     delete[] temp;
//   }

//   cout << "My rank is " << rank << ", and my program name is " 
//        << prog_name << " which is " << len << " chars long." << endl;
//   cout << "I'm going to load data from the file '" << inputfile << "'."
//        << endl;

  cout << "phread version " << PACKAGE_VERSION << " starting. " 
       << "Rank " << rank << " of " << size << "." << endl;

  // Parse the input script (each process will just load it's own file
  // for now. ) 

  // Subdomain the grid among the available processors and have each
  // processor set up its grid.

  SimpleSDAlg dd; // Domain decomposition algorithm. Maybe it should
                  // be static? 
  

  // THIS STUFF HERE IS TEMPORARY
  Grid grid; 
  GridInfo info_g;
  
  info_g.global_dimx_ = info_g.dimx_ = 100;
  info_g.global_dimy_ = info_g.dimy_ = 100;
  info_g.global_dimz_ = info_g.dimz_ = 100;
  info_g.deltax_ = info_g.deltay_ = info_g.deltaz_ = 18.75e-9;
  info_g.deltat_ = 36e-18;
  info_g.start_x_ = info_g.start_y_ = info_g.start_z_ = 0;

  info_g.set_boundary(FRONT, EWALL);
  info_g.set_boundary(BACK, EWALL);
  info_g.set_boundary(LEFT, EWALL);
  info_g.set_boundary(RIGHT, EWALL);
  info_g.set_boundary(BOTTOM, EWALL);
  info_g.set_boundary(TOP, EWALL);
  
  GridInfo info = dd.decompose_domain(rank, size, info_g);

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

  // Excitation
  Gaussm ex;
  ex.set_parameters(1, 500e12, 300e12);
  ex.set_region(50, 50, 50, 50, 50, 50);
  ex.set_polarization(0.0, 1.0, 0.0);

  // Results
  point_t p;
  p.x = 45;
  p.y = 50;
  p.z = 50;
  PointResult res1;
  PointResult res2;
  PointResult res3;
  res1.set_point(p);
  p.x = 50;
  res2.set_point(p);
  p.x = 55;
  res3.set_point(p);

  AsciiDataWriter adw1(rank, size);
  adw1.set_filename("t_ey_45_2.txt");
  adw1.add_variable(res1);
  
  AsciiDataWriter adw2(rank, size);
  adw2.set_filename("t_ey_50_2.txt");
  adw2.add_variable(res2);

  AsciiDataWriter adw3(rank, size);
  adw3.set_filename("t_ey_55_2.txt");
  adw3.add_variable(res3);

  grid.set_define_mode(false);
  
  // Main loop
  unsigned int num_time_steps = 7;
  unsigned int ts = 0;

  //ex.excite(grid, ts, BOTH);
  grid.apply_boundaries();
    
  cout << "main loop begins." << endl;  
  for (ts = 1; ts < num_time_steps; ts++) {
    cout << "phred, time step " << ts << ", excitation: " 
         << ex.source_function(grid, ts) << endl;

    // Fields update
    grid.update_fields();

    // Total / Scattered excitation

    // Excitations
    ex.excite(grid, ts, E);
    ex.excite(grid, ts, H);

    // Boundary condition application
    grid.apply_boundaries();

    // Total / Scattered field interface confditions

    // Results
    adw1.handle_data(ts, res1.get_result(grid, ts));
    adw2.handle_data(ts, res2.get_result(grid, ts));
    adw3.handle_data(ts, res3.get_result(grid, ts));
  }

  cout << "phred is phinished." << endl;

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
");
#else
  printf ("Usage: %s FILENAME\n", program_name);
  printf("\
FILENAME is the file containing the description of the proble to simulate.\n\
");
#endif

  exit (status);
}
