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
#include "NetCDFDataWriter.hh"
#include "PlaneResult.hh"
#include "PointDFTResult.hh"
#include "SourceDFTResult.hh"
#include "SourceTimeResult.hh"
#include "Excitation.hh"
#include "BartlettExcitation.hh"
#include "FDTD.hh"
#include "Box.hh"

#ifdef USE_PY_BINDINGS
#include <Python.h>
#include "python-bindings/PyInterpreter.hh"
#endif

static void usage (int status);

#ifdef HAVE_LIBPOPT

/* popt plays way nicer with MPI than getopt. Trust me. */
#include <popt.h>

static struct poptOption options[] = 
  {
    {"help", 'h', POPT_ARG_NONE, 0, 'h'},
    {"version", 'V', POPT_ARG_NONE, 0, 'V'},
    {"interactive", 'i', POPT_ARG_NONE, 0, 'i'},
    {"file", 'f', POPT_ARG_STRING, 0, 'f'},
    {0, 0, 0, 0, 0}
  };
#endif

static int decode_switches (int argc, char **argv);

// Ugly globals
string inputfile;
const char *program_name;
bool interactive;

// Test runs
static void point_test(int rank, int size);
static void pml_test(int rank, int size);

// MAIN!
int
main (int argc, char **argv)
{
  int i, rank, size, len;
  string prog_name;
  char *temp;
  interactive = false;

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

  cout << "phread version " << PACKAGE_VERSION << " starting on " 
       << "rank " << rank << " of " << size << " processes." << endl;

  // Parse the input script (each process will just load it's own file
  // for now. ) 

  // Subdomain the grid among the available processors and have each
  // processor set up its grid.

  try {
#ifdef USE_PY_BINDINGS
    if (interactive) {
      PyInterpreter interp().run(rank, size);
    } 
#endif
    if (!interactive) {
      // If the file extension is .py, run it in the python interpreter. 

      // If there is no extension, look for files named fn.in, fn.str,
      // and fn.out, per Jan's grammer, and load them instead. 

      // TESTS, TEMPORARY
      //point_test(rank, size);
      pml_test(rank, size);
    }
  } catch (const std::exception &e) {
    cout << "Caught exception: " << e.what() << endl;
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

    case 'i':
      interactive = true;
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

  
  return 0;
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
  -f, --filename             filename to read problem description from\n");
#ifdef USE_PY_BINDINGS
  printf("  -i, --interactive          start an interactive Python interpreter on\n\
                             rank zero if that process is attached to a \n\
                             terminal. Commands will be mirrored to\n\
                             interpreters running on the other ranks.");
#endif
  printf("  -h, --help                 display this help and exit\n\
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


// Test runs
static void point_test(int rank, int size)
{
  FDTD fdtd;
  
  fdtd.set_grid_size(100, 100, 100);
  fdtd.set_grid_deltas(18.75e-9, 18.75e-9, 18.75e-9);

  Ewall ewall;
  fdtd.set_boundary(FRONT, &ewall);
  fdtd.set_boundary(BACK, &ewall);
  fdtd.set_boundary(BOTTOM, &ewall);
  fdtd.set_boundary(TOP, &ewall);
  fdtd.set_boundary(LEFT, &ewall);
  fdtd.set_boundary(RIGHT, &ewall);
  
  MaterialLib mats; 
  Material mat; // defaults to free space
  mats.add_material(mat);
  
  // The library stores copies. 
  mat.set_epsilon(2.2);
  mat.set_name("Substrate");
  mats.add_material(mat);

  fdtd.load_materials(mats);

  // Global coordinates. 
  Box all;
  all.set_region(0, 100, 0, 100, 0, 100);
  all.set_material_id(1);

  Box diel;
  diel.set_region(40, 60, 40, 60, 40, 60);
  diel.set_material_id(2);

  fdtd.add_geometry(&all);
  fdtd.add_geometry(&diel);

  // Excitation
  Gaussm gm;
  gm.set_parameters(1, 500e12, 300e12);

  BartlettExcitation ex(&gm);
  ex.set_soft(false);
  ex.set_region(50, 50, 50, 50, 50, 50);
  ex.set_polarization(0.0, 1.0, 0.0);
  
  fdtd.add_excitation("modgauss", &ex);

  // Results
  point_t p;
  p.x = 50;
  p.y = 50;
  p.z = 50;
  PointResult res1(p);
  PointDFTResult pdft(100e12, 600e12, 50);
  pdft.set_point(p);

  fdtd.add_result("res1", &res1);
  fdtd.add_result("pdft", &pdft);

  AsciiDataWriter adw1(rank, size);
  adw1.set_filename("t_field_50.txt");

  AsciiDataWriter adw2(rank, size);
  adw2.set_filename("t_field_dft_50.txt");

  fdtd.add_datawriter("adw1", &adw1);
  fdtd.add_datawriter("adw2", &adw2);
  fdtd.map_result_to_datawriter("res1", "adw1");
  fdtd.map_result_to_datawriter("pdft", "adw2");

  NetCDFDataWriter ncdw(rank, size);
  ncdw.set_filename("yz_plane.nc");

  fdtd.add_datawriter("ncdw", &ncdw);

  PlaneResult pr1;
  pr1.set_name("yzplane");
  pr1.set_plane(p, BACK);

  fdtd.add_result("pr1", &pr1);
  fdtd.map_result_to_datawriter("pr1", "ncdw");

  SourceDFTResult sdftr(gm, 100e12, 600e12, 50);
  sdftr.set_time_param(0, 100, 0);
  fdtd.add_result("sdftr", &sdftr);

  AsciiDataWriter adw5(rank, size);
  adw5.set_filename("src_dft.txt");
  fdtd.add_datawriter("adw5", &adw5);
  
  fdtd.map_result_to_datawriter("sdftr", "adw5");

  fdtd.run(rank, size, 100);
}

static void pml_test(int rank, int size)
{
  FDTD fdtd;
  
  fdtd.set_grid_size(100, 50, 60);
  fdtd.set_grid_deltas(18.75e-9, 18.75e-9, 18.75e-9);

  Pml front(VP, 1.0), back(VP, 1.0), left(VP, 1.0), right(VP, 1.0),
    top(VP, 1.0), bottom(VP, 1.0);
  front.set_thickness(4);
  back.set_thickness(4);
  left.set_thickness(4);
  right.set_thickness(4);
  top.set_thickness(4);
  bottom.set_thickness(4);

  fdtd.set_boundary(FRONT, &front);
  fdtd.set_boundary(BACK, &back);
  fdtd.set_boundary(BOTTOM, &bottom);
  fdtd.set_boundary(TOP, &top);
  fdtd.set_boundary(LEFT, &left);
  fdtd.set_boundary(RIGHT, &right);

  MaterialLib mats; 
  Material mat; // defaults to free space
  mats.add_material(mat);
  
  // The library stores copies. 
  mat.set_epsilon(2.2);
  mat.set_name("Substrate");
  mats.add_material(mat);

  fdtd.load_materials(mats);

  // Global coordinates. 
  Box all;
  all.set_region(0, 100, 0, 50, 0, 60);
  all.set_material_id(1);

  Box metal1;
  metal1.set_region(45, 55, 5, 46, 0, 25);
  metal1.set_material_id(0);

  Box metal2;
  metal2.set_region(45, 55, 5, 46, 35, 60);
  metal2.set_material_id(0);

  fdtd.add_geometry(&all);
  fdtd.add_geometry(&metal1);
  fdtd.add_geometry(&metal2);

  // Excitation
  Gaussm gm;
  gm.set_parameters(1, 500e12, 300e12);

  Excitation ex(&gm);
  //BartlettExcitation ex(gm);
  ex.set_soft(true);
  ex.set_region(15, 15, 10, 40, 20, 40);
  ex.set_polarization(0.0, 1.0, 0.0);

  fdtd.add_excitation("modgauss", &ex);

  // Results
  point_t p;
  p.x = 20;
  p.y = 25;
  p.z = 30;
  PointResult res1(p);
  PointDFTResult pdft(100e12, 600e12, 50);
  pdft.set_point(p);

  fdtd.add_result("res1", &res1);
  fdtd.add_result("pdft", &pdft);
  
  AsciiDataWriter adw1(rank, size);
  adw1.set_filename("t_field_20.txt");

  AsciiDataWriter adw2(rank, size);
  adw2.set_filename("t_field_dft_20.txt");

  fdtd.add_datawriter("adw1", &adw1);
  fdtd.add_datawriter("adw2", &adw2);
  fdtd.map_result_to_datawriter("res1", "adw1");
  fdtd.map_result_to_datawriter("pdft", "adw2");

  NetCDFDataWriter ncdw(rank, size);
  ncdw.set_filename("yz_plane.nc");

  fdtd.add_datawriter("ncdw", &ncdw);

  PlaneResult pr1;
  pr1.set_name("hx-xzplane");
  pr1.set_plane(p, LEFT);
  pr1.set_field(FC_HX);
  
  PlaneResult pr2;
  pr2.set_name("ey-xzplane");
  pr2.set_plane(p, LEFT);
  pr2.set_field(FC_EY);

   PlaneResult pr3;
   pr3.set_name("hz-xzplane");
   pr3.set_plane(p, LEFT);
   pr3.set_field(FC_HZ);

   PlaneResult pr4;
   pr4.set_name("hx-yzplane4");
   pr4.set_plane(p, FRONT);
   pr4.set_field(FC_HX);

   PlaneResult pr5;
   pr5.set_name("ey-xyplane5");
   pr5.set_plane(p, BOTTOM);
   pr5.set_field(FC_EY);

   fdtd.add_result("pr1", &pr1);
   fdtd.add_result("pr2", &pr2);
   fdtd.add_result("pr3", &pr3);
   fdtd.add_result("pr4", &pr4);
   fdtd.add_result("pr5", &pr5);

   fdtd.map_result_to_datawriter("pr1", "ncdw");
   fdtd.map_result_to_datawriter("pr2", "ncdw");
   fdtd.map_result_to_datawriter("pr3", "ncdw");
   fdtd.map_result_to_datawriter("pr4", "ncdw");
   fdtd.map_result_to_datawriter("pr5", "ncdw");

   SourceDFTResult sdftr(gm, 100e12, 600e12, 50);
   sdftr.set_time_param(0, 500, 0);
   fdtd.add_result("sdftr", &sdftr);

   AsciiDataWriter adw5(rank, size);
   adw5.set_filename("src_dft.txt");
   fdtd.add_datawriter("adw5", &adw5);

   fdtd.map_result_to_datawriter("sdftr", "adw5");

   SourceTimeResult srctr(gm);
   AsciiDataWriter adw8(rank, size);
   adw8.set_filename("src.txt");

   fdtd.add_result("srctr", &srctr);
   fdtd.add_datawriter("adw8", &adw8);

   fdtd.map_result_to_datawriter("srctr", "adw8");

   fdtd.run(rank, size, 100);
}

