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

#ifdef USE_PY_BINDINGS
#include <Python.h>
#endif

static void usage (int status);

#ifdef HAVE_LIBPOPT

/* popt plays way nicer with MPI than getopt. Trust me. */
#include <popt.h>

static struct poptOption options[] = 
  {
    {"help", 'h', POPT_ARG_NONE, 0, 'h'},
    {"version", 'V', POPT_ARG_NONE, 0, 'V'},
    {"verbose", 'v', POPT_ARG_NONE, 0, 'v'},
    {"file", 'f', POPT_ARG_STRING, 0, 'f'},
    {0, 0, 0, 0, 0}
  };
#endif

static int decode_switches (int argc, char **argv);

// Ugly globals
string inputfile;
const char *program_name;

// Some test declarations
static void test_yz_plane(field_t *ptr, MPI_Datatype t, int xpos, int rank);
static void test2(field_t *ptr, MPI_Datatype t);
static void test3_yz(Grid &grid, int i);
static void test4_z_vec(field_t *ptr, MPI_Datatype t);
static void test5_y_vec(field_t *ptr, MPI_Datatype t);
static void test6_x_vec(field_t *ptr, MPI_Datatype t);

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

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

#ifdef USE_PY_BINDINGS
  Py_Initialize();
#endif 

  if (rank == 0)
  {
    prog_name = argv[0];
    len = prog_name.size();
    program_name = prog_name.c_str();
    
    // MPI implementations are not required to distribute command line
    // args, although MPICH does.
    i = decode_switches (argc, argv);
  } 

  cout << "phread version " << PACKAGE_VERSION << " starting. " 
       << "Rank " << rank << " of " << size << "." << endl;

  // Parse the input script (each process will just load it's own file
  // for now. ) 

  // Subdomain the grid among the available processors and have each
  // processor set up its grid.

  try {
    //point_test(rank, size);
    pml_test(rank, size);
  } catch (const std::exception &e) {
    cout << "Caught exception: " << e.what() << endl;
  }

  cout << "phred is phinished." << endl;

#ifdef USE_PY_BINDINGS
  Py_Finalize();
#endif 

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


// Test runs
static void point_test(int rank, int size)
{
  SimpleSDAlg dd; // Domain decomposition algorithm. Maybe it should
                  // be static? 
  

  // THIS STUFF HERE IS TEMPORARY
  Grid grid; 
  GridInfo info_g;
  
  info_g.global_dimx_ = info_g.dimx_ = 100;
  info_g.global_dimy_ = info_g.dimy_ = 100;
  info_g.global_dimz_ = info_g.dimz_ = 100;
  info_g.deltax_ = 18.75e-9; 
  info_g.deltay_ = 18.75e-9; 
  info_g.deltaz_ = 18.75e-9;
  info_g.deltat_ = 36e-18;
  info_g.start_x_ = info_g.start_y_ = info_g.start_z_ = 0;


  info_g.set_boundary(FRONT, EWALL);
  info_g.set_boundary(BACK, EWALL);
  info_g.set_boundary(BOTTOM, EWALL);
  info_g.set_boundary(TOP, EWALL);
  info_g.set_boundary(LEFT, EWALL);
  info_g.set_boundary(RIGHT, EWALL);

  
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
  ex.set_soft(false);
  ex.set_parameters(1, 500e12, 300e12);
  ex.set_region(50, 50, 50, 50, 50, 50);
  ex.set_polarization(0.0, 1.0, 0.0);

  // Results
  point_t p;
  p.x = 50;
  p.y = 50;
  p.z = 50;
  PointResult res1;
  res1.set_point(p);

  AsciiDataWriter adw1(rank, size);
  adw1.set_filename("t_field_50.txt");
  adw1.add_variable(res1);

  PointDFTResult pdft(100e12, 600e12, 50);
  pdft.set_point(p);

  AsciiDataWriter adw2(rank, size);
  adw2.set_filename("t_field_dft_50.txt");
  adw2.add_variable(pdft);
  
  //AsciiDataWriter adw4(rank, size);
  NetCDFDataWriter ncdw(rank, size);
  ncdw.set_filename("yz_plane.nc");
  ncdw.init();
  PlaneResult pr1;
  pr1.set_name("yzplane");
  pr1.set_plane(p, BACK);
  pr1.set_size(grid.get_ldy(), grid.get_ldz());
  //adw4.set_filename("yz_plane.txt");
  //adw4.add_variable(pr1);
  ncdw.add_variable(pr1);

  SourceDFTResult sdftr(ex, 100e12, 600e12, 50);
  sdftr.set_time_param(0, 100, 0);

  AsciiDataWriter adw5(rank, size);
  adw5.set_filename("src_dft.txt");
  adw5.add_variable(sdftr);

  grid.set_define_mode(false);
  
  // Main loop
  unsigned int num_time_steps = 101;
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

    //test6_x_vec(grid.get_face_start(BACK, EY, 0), grid.get_x_vector_dt());
    //test5_y_vec(grid.get_face_start(BACK, EY, 0), grid.get_y_vector_dt());
    //test4_z_vec(grid.get_face_start(BACK, EY, 0), grid.get_z_vector_dt());
    //test(grid.get_face_start(FRONT, EY, p2), grid.get_plane_dt(FRONT));
    //cout << "\nTest2: " << endl;
    //test2(grid.get_face_start(BACK, EY, 2), grid.get_plane_dt(FRONT));
    //test3_yz(grid, 2);
    //cout << "\nTest trans: " << endl;
    //test_trans(grid.get_face_start(BACK, EY), grid.get_plane_dt(BACK));
    //if (rank == 0)
    //  test3(grid, 25);
    //else
    //test3(grid, 0);

    //if (rank == 0)
    //test_yz_plane(grid.get_face_start(BACK, EY, 22), 
    //                grid.get_plane_dt(BACK), 0, rank);
    //else if (rank == 1)
    //  test_yz_plane(grid.get_face_start(BACK, EY, 0), 
    //                grid.get_plane_dt(BACK), 0, rank);

    // Boundary condition application
    grid.apply_boundaries();

    // Total / Scattered field interface confditions

    // Results

    adw1.handle_data(ts, res1.get_result(grid, ts));
    adw2.handle_data(ts, pdft.get_result(grid, ts));
    //adw4.handle_data(ts, pr1.get_result(grid, ts));
    adw5.handle_data(ts, sdftr.get_result(grid, ts));
    ncdw.handle_data(ts, pr1.get_result(grid, ts));
  }
}

static void pml_test(int rank, int size)
{
  SimpleSDAlg dd; // Domain decomposition algorithm. Maybe it should
                  // be static? 
  

  // THIS STUFF HERE IS TEMPORARY
  Grid grid; 
  GridInfo info_g;

  info_g.global_dimx_ = info_g.dimx_ = 40;
  info_g.global_dimy_ = info_g.dimy_ = 50;
  info_g.global_dimz_ = info_g.dimz_ = 60;
  info_g.deltax_ = 18.75e-9; 
  info_g.deltay_ = 18.75e-9; 
  info_g.deltaz_ = 18.75e-9;
  info_g.deltat_ = 36e-18;
  info_g.start_x_ = info_g.start_y_ = info_g.start_z_ = 0;

  Pml *pml = dynamic_cast<Pml *>(info_g.set_boundary(FRONT, PML));
  pml->set_thickness(4);
  pml->set_variation(VP);
  pml->set_nrml_refl(1.0);

  pml = dynamic_cast<Pml *>(info_g.set_boundary(BACK, PML));
  pml->set_thickness(4);
  pml->set_variation(VP);
  pml->set_nrml_refl(1.0);

  pml = dynamic_cast<Pml *>(info_g.set_boundary(LEFT, PML));
  pml->set_thickness(4);
  pml->set_variation(VP);
  pml->set_nrml_refl(1.0);

  pml = dynamic_cast<Pml *>(info_g.set_boundary(RIGHT, PML));
  pml->set_thickness(4);
  pml->set_variation(VP);
  pml->set_nrml_refl(1.0);

  pml = dynamic_cast<Pml *>(info_g.set_boundary(BOTTOM, PML));
  pml->set_thickness(4);
  pml->set_variation(VP);
  pml->set_nrml_refl(1.0);

  pml = dynamic_cast<Pml *>(info_g.set_boundary(TOP, PML));
  pml->set_thickness(4);
  pml->set_variation(VP);
  pml->set_nrml_refl(1.0);

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
  grid.define_box(0, 40, 0, 50, 0, 60, 1);

  // Excitation
  Gaussm ex;
  ex.set_soft(false);
  ex.set_parameters(1, 500e12, 300e12);
  ex.set_region(20, 20, 25, 25, 30, 30);
  ex.set_polarization(0.0, 1.0, 0.0);

  // Results
  point_t p;
  p.x = 20;
  p.y = 25;
  p.z = 30;
  PointResult res1;
  res1.set_point(p);
  
  AsciiDataWriter adw1(rank, size);
  adw1.set_filename("t_field_20.txt");
  adw1.add_variable(res1);

  PointDFTResult pdft(100e12, 600e12, 50);
  pdft.set_point(p);

  AsciiDataWriter adw2(rank, size);
  adw2.set_filename("t_field_dft_20.txt");
  adw2.add_variable(pdft);

  point_t p2;
  p2.x = 5;
  p2.y = 25;
  p2.z = 30;
  PointResult res2;
  res2.set_point(p2);
  
  AsciiDataWriter adw3(rank, size);
  adw3.set_filename("t_field_5.txt");
  adw3.add_variable(res2);

  p2.x = 10;
  PointResult res3;
  res3.set_point(p2);

  AsciiDataWriter adw4(rank, size);
  adw4.set_filename("t_field_10.txt");
  adw4.add_variable(res3);

  p2.x = 11;
  PointResult res4;
  res4.set_point(p2);
  
  AsciiDataWriter adw6(rank, size);
  adw6.set_filename("t_field_11.txt");
  adw6.add_variable(res4);

  p2.x = 9;
  PointResult res5;
  res5.set_point(p2);
  
  AsciiDataWriter adw7(rank, size);
  adw7.set_filename("t_field_9.txt");
  adw7.add_variable(res5);


//   //AsciiDataWriter adw4(rank, size);
   NetCDFDataWriter ncdw(rank, size);
   ncdw.set_filename("yz_plane.nc");
   ncdw.init();
   PlaneResult pr1;
   pr1.set_name("yzplane");
   pr1.set_plane(p, BACK);
   pr1.set_size(grid.get_ldy(), grid.get_ldz());

   PlaneResult pr2;
   pr2.set_name("xzplane");
   pr2.set_plane(p, LEFT);
   pr2.set_size(grid.get_ldx(), grid.get_ldz());

   PlaneResult pr3;
   pr3.set_name("xyplane");
   pr3.set_plane(p, BOTTOM);
   pr3.set_size(grid.get_ldz(), grid.get_ldy());

   //adw4.set_filename("yz_plane.txt");
   //adw4.add_variable(pr1);
   ncdw.add_variable(pr1);
   ncdw.add_variable(pr2);
   ncdw.add_variable(pr3);

   SourceDFTResult sdftr(ex, 100e12, 600e12, 50);
   sdftr.set_time_param(0, 500, 0);

   AsciiDataWriter adw5(rank, size);
   adw5.set_filename("src_dft.txt");
   adw5.add_variable(sdftr);

   SourceTimeResult srctr(ex);
   AsciiDataWriter adw8(rank, size);
   adw8.set_filename("src.txt");
   adw8.add_variable(srctr);
   

  grid.set_define_mode(false);
  
  // Main loop
  unsigned int num_time_steps = 101;
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

    //test6_x_vec(grid.get_face_start(BACK, EY, 0), grid.get_x_vector_dt());
    //test5_y_vec(grid.get_face_start(BACK, EY, 0), grid.get_y_vector_dt());
    //test4_z_vec(grid.get_face_start(BACK, EY, 0), grid.get_z_vector_dt());
    //test(grid.get_face_start(FRONT, EY, p2), grid.get_plane_dt(FRONT));
    //cout << "\nTest2: " << endl;
    //test2(grid.get_face_start(BACK, EY, 2), grid.get_plane_dt(FRONT));
    //test3_yz(grid, 2);
    //cout << "\nTest trans: " << endl;
    //test_trans(grid.get_face_start(BACK, EY), grid.get_plane_dt(BACK));
    //if (rank == 0)
    //  test3(grid, 25);
    //else
    //test3(grid, 0);

    //if (rank == 0)
    //test_yz_plane(grid.get_face_start(BACK, EY, 22), 
    //                grid.get_plane_dt(BACK), 0, rank);
    //else if (rank == 1)
    //  test_yz_plane(grid.get_face_start(BACK, EY, 0), 
    //                grid.get_plane_dt(BACK), 0, rank);

    // Boundary condition application
    grid.apply_boundaries();

    // Total / Scattered field interface confditions

    // Results

     adw1.handle_data(ts, res1.get_result(grid, ts));
     adw3.handle_data(ts, res2.get_result(grid, ts));
     adw4.handle_data(ts, res3.get_result(grid, ts));
     adw6.handle_data(ts, res4.get_result(grid, ts));
     adw7.handle_data(ts, res5.get_result(grid, ts));
     adw8.handle_data(ts, srctr.get_result(grid, ts));
     adw2.handle_data(ts, pdft.get_result(grid, ts));
     //adw4.handle_data(ts, pr1.get_result(grid, ts));
     adw5.handle_data(ts, sdftr.get_result(grid, ts));
     ncdw.handle_data(ts, pr1.get_result(grid, ts));
     ncdw.handle_data(ts, pr2.get_result(grid, ts));
     //ncdw.handle_data(ts, pr3.get_result(grid, ts));


  }
}


// TESTS
static void test_yz_plane(field_t *ptr, MPI_Datatype t, int xpos, int rank)
{
  int sz, pos, floats,j,k,i;
  field_t *r, *sum, *sum2, **sums;
  MPI_Status status;
  MPI_Datatype rt;

  MPI_Type_size(t, &sz);

  floats = sz / sizeof(field_t);

  r = new field_t[floats];
  
  cout << "Testing MPI yz plane datatype at x = " << xpos << endl;

  ptr += 0 + (0 + xpos*50)*50;

  MPI_Type_contiguous(floats, MPI_FLOAT, &rt);
  MPI_Type_commit(&rt);

  MPI_Sendrecv(static_cast<void *>(ptr), 1, t, rank, 1,
               static_cast<void *>(r), 1, rt, rank, 1,
               MPI_COMM_WORLD, &status);  
  
  int idx = 0;
  for (i = 0; i < 50; i++) {
    for (j = 0; j < 50; j++) {
      cout.width(4);
      cout << r[idx++];
    }
    cout << endl;
  }

  cout << endl;

  delete[] r;    
  
}

// This one works
static void test2(field_t *ptr, MPI_Datatype t)
{
  int sz, pos, floats;
  field_t *r;

  MPI_Type_size(t, &sz);

  floats = sz / sizeof(field_t);

  cout << "Simple loop: Subdomain plane is some " << sz << " bytes in size, or "
       << floats << " floats..." << endl;

  unsigned int i, j, k;
  for (j=0; j < 50; j++) {
    for (k=0; k < 50; k++) {
      cout << ptr[k + j*50] << " ";
      if ((k+j*50) % 50 == 0) 
        cout << endl;
    }
  }
  
}

// This one works
static void test3_yz(Grid &grid, int i)
{
  int sz, pos, floats;

  unsigned int j, k;

  cout << "Another simple loop, calling grid functions..." << endl;
  for (j=0; j < 50; j++) {
    for (k=0; k < 50; k++) {
      cout << grid.get_ey(i, j, k) << " ";
      if ((k+j*50) % 50 == 0) 
        cout << endl;
    }
  }
}

// This one works
static void test4_z_vec(field_t *ptr, MPI_Datatype t)
{
  int sz, pos, floats, k, j;
  field_t *r;

  MPI_Type_size(t, &sz);

  floats = sz / sizeof(field_t);

  cout << "Z: Vector is some " << sz << " bytes in size, or "
       << floats << " floats..." << endl;

  ptr += (25 + 4*50)*50;

  for (k=0; k < 50; k++) {
    cout << ptr[k] << " ";
    if ((k % 50 == 0))
      cout << endl;
  }
  cout << endl << "Z: And the MPI way: " << endl;

  r = new field_t[floats];
  pos = 0;
  MPI_Pack(static_cast<void *>(ptr), 1, t, 
           static_cast<void *>(r), sz, &pos, MPI_COMM_WORLD);

  cout << "Z: And here it is, pos = " << pos << "... " << endl;

  for (int i = 0; i < floats; i++) {
    cout << r[i] << " ";
    if (i % 50 == 0)
      cout << endl;
  }
  cout << endl;

  delete[] r;
}

// This one works
static void test5_y_vec(field_t *ptr, MPI_Datatype t)
{
  int sz, pos, floats, k, j;
  field_t *r;

  MPI_Type_size(t, &sz);

  floats = sz / sizeof(field_t);

  cout << "Y: Vector is some " << sz << " bytes in size, or "
       << floats << " floats..." << endl;

  ptr += 25 + (0 + 4*50)*50;

  for (j=0; j < 50; j++) {
    cout << ptr[j*50] << " ";
    if ((j % 50 == 0))
      cout << endl;
  }
  cout << endl;

  cout << endl << "And the MPI way: " << endl;

  r = new field_t[floats];
  pos = 0;
  MPI_Pack(static_cast<void *>(ptr), 1, t, 
           static_cast<void *>(r), sz, &pos, MPI_COMM_WORLD);

  cout << "Y: And here it is, pos = " << pos << "... " << endl;

  for (int i = 0; i < floats; i++) {
    cout << r[i] << " ";
    if (i % 50 == 0)
      cout << endl;
  }
  cout << endl;

  delete[] r;
}

// This one works
static void test6_x_vec(field_t *ptr, MPI_Datatype t)
{
  int sz, pos, floats, k, j, i;
  field_t *r;

  MPI_Type_size(t, &sz);

  floats = sz / sizeof(field_t);

  cout << "X: Vector is some " << sz << " bytes in size, or "
       << floats << " floats..." << endl;

  ptr += 25 + (25 + 0*50)*50;

  for (i=0; i < 50; i++) {
    cout << ptr[i*50*50] << " ";
    if ((i % 50 == 0))
      cout << endl;
  }
  cout << endl;

  cout << endl << "And the MPI way: " << endl;

  r = new field_t[floats];
  pos = 0;
  MPI_Pack(static_cast<void *>(ptr), 1, t, 
           static_cast<void *>(r), sz, &pos, MPI_COMM_WORLD);

  cout << "X: And here it is, pos = " << pos << "... " << endl;

  for (int i = 0; i < floats; i++) {
    cout << r[i] << " ";
    if (i % 50 == 0)
      cout << endl;
  }
  cout << endl;

  delete[] r;
}
