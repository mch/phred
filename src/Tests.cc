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

// In order to actually do stuff...
#include "MaterialLib.hh"
#include "Grid.hh"
#include "SimpleSDAlg.hh"
#include "Signals/Gaussm.hh"
#include "Boundaries/Ewall.hh"
#include "Boundaries/Hwall.hh"
#include "Boundaries/UPml.hh"
#include "Results/PointResult.hh"
#include "Results/GridResult.hh"
#include "Results/FarfieldResult.hh"
#include "DataWriters/AsciiDataWriter.hh"
#include "DataWriters/MatlabDataWriter.hh"
#include "DataWriters/NetCDFDataWriter.hh"
#include "Results/PlaneResult.hh"
#include "Results/PointDFTResult.hh"
#include "Results/PowerResult.hh"
#include "Results/SignalDFTResult.hh"
#include "Results/SignalTimeResult.hh"
#include "Excitations/Excitation.hh"
#include "Excitations/BartlettExcitation.hh"
#include "Excitations/GaussWindExcitation.hh"
#include "Excitations/WaveguideExcitation.hh"
#include "CSG/CSGCylinder.hh"
#include "CSG/CSGDifference.hh"
#include "FDTD.hh"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <string>
#include <iostream>

using namespace std;

/**
 * Million node benchmark
 */ 
void mn_benchmark()
{
  FDTD fdtd;
  fdtd.set_grid_deltas(1e-2, 1e-2, 1e-2);
  fdtd.set_grid_size(1,1,1);

  fdtd.set_time_steps(100);

#ifdef USE_OPENMP
   // Test the OpenMP
   time_t start, now;
   clock_t cpu_start, cpu_now;

   for (int numthreads = 1; numthreads <= omp_get_max_threads(); numthreads++)
     {
       omp_set_num_threads(numthreads);
       
       cout << "100x100x100 million node 100 time step benchmark on "
            << numthreads << " of " << omp_get_max_threads() << "...\n";

       start = time(NULL);
       cpu_start = clock();
       fdtd.run();
       now = time(NULL);
       cpu_now = clock();
       
       cout << numthreads << " of " 
	    << omp_get_max_threads() << " threads took " 
	    << now - start << " wall clock seconds, and "
	    << (cpu_now - cpu_start) / static_cast<double>(CLOCKS_PER_SEC)
	    << " cpu seconds." << endl;
     }

#else
   cout << "100x100x100 million node 100 time step benchmark\n";
   fdtd.run();
#endif

}

/**
 * Variable size benchmark
 */ 
void var_benchmark(unsigned int x_cells, unsigned int y_cells, 
                   unsigned int z_cells)
{
  if (x_cells < 100 || y_cells < 100 || z_cells < 100)
  {
    cout << "Don't waste your time. Use a bigger problem for benchmarking.\n";
    return;
  }

  FDTD fdtd;
  fdtd.set_grid_deltas(1.0 / static_cast<float>(x_cells), 
                       1.0 / static_cast<float>(y_cells),
                       1.0 / static_cast<float>(z_cells));
  fdtd.set_grid_size(1,1,1);

  fdtd.set_time_steps(100);

#ifdef USE_OPENMP
   // Test the OpenMP
   time_t start, now;
   clock_t cpu_start, cpu_now;

   for (int numthreads = 1; numthreads <= omp_get_max_threads(); numthreads++)
     {
       omp_set_num_threads(numthreads);
       
       cout << "Variable size, 100 time step benchmark on "
            << numthreads << " of " << omp_get_max_threads() << "...\n";

       start = time(NULL);
       cpu_start = clock();
       fdtd.run();
       now = time(NULL);
       cpu_now = clock();
       
       cout << numthreads << " of " 
	    << omp_get_max_threads() << " threads took " 
	    << now - start << " wall clock seconds, and "
	    << (cpu_now - cpu_start) / static_cast<double>(CLOCKS_PER_SEC)
	    << " cpu seconds." << endl;
     }

#else
   cout << "Variable size, 100 time step benchmark\n";
   fdtd.run();
#endif
 
}

/**
 * Simulates a single hole in a plate of perfect conductor
 */ 
void hole()
{
  cout << "Simulating a single hole in a plate of perfect conductor..." 
       << endl;

  string prefix = "hole_";

  // Time steps is recalculated based on the length of the excitation.
  unsigned int time_steps = 3000;

  float deltax = 5e-9;
  float deltay = 5e-9;
  float deltaz = 5e-9;

  float gridx = 800e-9;
  float gridy = 800e-9;
  float gridz = 800e-9;

  float plate_thickness = 100e-9;
  float hole_radius = 100e-9;

  // Excitation parameters
  float ex_ampl = 1.0;
  float ex_freq_size = 100e12;
  float ex_centre_f = 500e12;
  float ex_offset = -250e-9;

  // DFT Parameters
  float dft_low = 300e12;
  float dft_high = 750e12;
  unsigned int dft_num = 20;

  // GRID
  FDTD fdtd;

  fdtd.set_grid_deltas(deltax, deltay, deltaz);
  fdtd.set_grid_size(gridx, gridy, gridz);

  // MATERIAL
  shared_ptr<MaterialLib> mlib 
    = shared_ptr<MaterialLib>(new MaterialLib());

  Material mat;
  mat.set_epsilon(2.2);
  (*mlib).add_material("dielectric", mat);

  fdtd.load_materials(mlib);

  // BOUNDARIES
  for (int i = 0; i < 6; i++)
  {
    shared_ptr<Pml> bound = shared_ptr<Pml>(new Pml());
    (*bound).set_thickness(8);
    fdtd.set_boundary(static_cast<Face>(i), bound);
  }

  // EXCITATION
  shared_ptr<Gaussm> gm = shared_ptr<Gaussm>(new Gaussm());
  gm->set_parameters(ex_ampl, ex_freq_size, ex_centre_f);

  time_steps = static_cast<unsigned int>
    (ceil(gm->length() / fdtd.get_time_delta()));
  
  shared_ptr<GaussWindExcitation> ex 
    = shared_ptr<GaussWindExcitation>(new GaussWindExcitation(gm));
  shared_ptr<CSGBox> exbox 
    = shared_ptr<CSGBox>(new CSGBox());
  exbox->set_size(gridx - 16 * deltax, gridy - 16 * deltay, deltaz);
  exbox->set_centre(0, 0, ex_offset);
  ex->set_region(exbox);
  ex->set_soft(true);
  ex->set_type(E);
  ex->set_polarization(1,0,0);

  fdtd.add_excitation("fluffy", ex);

  // DATA WRITERS
//   shared_ptr<NetCDFDataWriter> ncdw 
//     = shared_ptr<NetCDFDataWriter>(new NetCDFDataWriter());

//   (*ncdw).set_filename(prefix + "planes.nc");
//   fdtd.add_datawriter("ncdw", ncdw);
  
  shared_ptr<MatlabDataWriter> mdw 
    = shared_ptr<MatlabDataWriter>(new MatlabDataWriter());

  (*mdw).set_filename(prefix + "power.mat");
  fdtd.add_datawriter("mdw", mdw);

  shared_ptr<GridResult> gridr
    = shared_ptr<GridResult>(new GridResult);

  // GRID RESULT: Should be disabled for the full problem
//   fdtd.add_result("grid", gridr);
//   fdtd.map_result_to_datawriter("grid", "ncdw");

  // INFORMATION ABOUT EXCIATION
  shared_ptr<SignalTimeResult> st
    = shared_ptr<SignalTimeResult>(new SignalTimeResult(*gm));
  fdtd.add_result("src", st);
  fdtd.map_result_to_datawriter("src", "mdw");

  shared_ptr<SignalDFTResult> sdft
    = shared_ptr<SignalDFTResult>(new SignalDFTResult(*gm, dft_low, 
                                                      dft_high, 
                                                      dft_num));
  fdtd.add_result("srcdft", sdft);
  fdtd.map_result_to_datawriter("srcdft", "mdw");

  // Incident side point
  shared_ptr<PointResult> pr1
    = shared_ptr<PointResult>(new PointResult());
  pr1->set_point(point(0, 0, -plate_thickness));
  fdtd.add_result("p1", pr1);

  shared_ptr<PointDFTResult> pr1dft 
    =  shared_ptr<PointDFTResult>(new PointDFTResult());
  pr1dft->set_point(point(0, 0, -plate_thickness));
  pr1dft->set_freq(dft_low, dft_high, dft_num);
  fdtd.add_result("p1dft", pr1dft);

  fdtd.map_result_to_datawriter("p1", "mdw");
  fdtd.map_result_to_datawriter("p1dft", "mdw");

  // Transmission side point
  shared_ptr<PointResult> pr2
    = shared_ptr<PointResult>(new PointResult());
  pr2->set_point(point(0, 0, plate_thickness));
  fdtd.add_result("p2", pr2);

  shared_ptr<PointDFTResult> pr2dft 
    =  shared_ptr<PointDFTResult>(new PointDFTResult());
  pr2dft->set_point(point(0, 0, plate_thickness));
  pr2dft->set_freq(dft_low, dft_high, dft_num);
  fdtd.add_result("p2dft", pr2dft);
  
  fdtd.map_result_to_datawriter("p2", "mdw");
  fdtd.map_result_to_datawriter("p2dft", "mdw");

  // Farfield measurements
  shared_ptr<CSGBox> ffbox
    = shared_ptr<CSGBox>(new CSGBox());
  ffbox->set_size(gridx, gridy, plate_thickness);
  
  // About the y axis
  shared_ptr<FarfieldResult> ffy
    = shared_ptr<FarfieldResult>(new FarfieldResult());
  ffy->set_freq(dft_low, dft_high, dft_num);
  ffy->set_region(ffbox);
  ffy->use_face(FRONT, false);
  ffy->use_face(BACK, false);
  ffy->use_face(LEFT, false);
  ffy->use_face(RIGHT, false);
  ffy->use_face(BOTTOM, false);
  ffy->use_face(TOP, true);  
  
  ffy->set_theta_degrees(-9, 9, 7);
  ffy->set_phi_degrees(0, 0, 1);

  // About the x axis
  shared_ptr<FarfieldResult> ffx
    = shared_ptr<FarfieldResult>(new FarfieldResult());
  ffx->set_freq(dft_low, dft_high, dft_num);
  ffx->set_region(ffbox);
  ffx->use_face(FRONT, false);
  ffx->use_face(BACK, false);
  ffx->use_face(LEFT, false);
  ffx->use_face(RIGHT, false);
  ffx->use_face(BOTTOM, false);
  ffx->use_face(TOP, true);  
  
  ffx->set_theta_degrees(-9, 9, 7);
  ffx->set_phi_degrees(90, 90, 1);

  fdtd.add_result("ffy", ffy);
  fdtd.add_result("ffx", ffx);
  
  fdtd.map_result_to_datawriter("ffy", "mdw");
  fdtd.map_result_to_datawriter("ffx", "mdw");
    
  shared_ptr<CSGBox> metal = shared_ptr<CSGBox>(new CSGBox());
  metal->set_size(gridx, gridy, plate_thickness);

  shared_ptr<CSGCylinder> hole 
    = shared_ptr<CSGCylinder>(new CSGCylinder());
  hole->set_radius(hole_radius);
  hole->set_height(2 * plate_thickness);
  
  shared_ptr<CSGDifference> plate 
    = shared_ptr<CSGDifference>(new CSGDifference(metal, hole));
  
  fdtd.add_object("PEC", plate);

  fdtd.set_time_steps(time_steps);
  fdtd.run();
}

void grooves_top()
{}

void grooves_bottom()
{}

void grooves_both()
{}

// Test runs
// void point_test(int rank, int size)
// {
//   FDTD fdtd;
  
//   fdtd.set_grid_size(100, 100, 400);
//   fdtd.set_grid_deltas(18.75e-9, 18.75e-9, 18.75e-9);

//   Ewall ewall;
//   fdtd.set_boundary(FRONT, &ewall);
//   fdtd.set_boundary(BACK, &ewall);
//   fdtd.set_boundary(BOTTOM, &ewall);
//   fdtd.set_boundary(TOP, &ewall);
//   fdtd.set_boundary(LEFT, &ewall);
//   fdtd.set_boundary(RIGHT, &ewall);
  
//   MaterialLib mats; 
//   Material mat; // defaults to free space
//   mats.add_material(mat);
  
//   // The library stores copies. 
//   mat.set_epsilon(2.2);
//   mat.set_name("Substrate");
//   mats.add_material(mat);

//   mat.set_epsilon(1);
//   mat.set_name("Silver");
//   //mat.set_collision_freq(57e12); // THz
//   //mat.set_plasma_freq(2 * PI * 2000e+12); // THz * 2 * pi
//   mat.set_collision_freq(1.4e+14);
//   mat.set_plasma_freq(2 * PI * 1.85e+15);
//   mats.add_material(mat);

//   fdtd.load_materials(mats);

//   // Global coordinates. 
//   Box all;
//   all.set_region(0, 100, 0, 100, 0, 100);
//   all.set_material_id(1);

//   Box diel;
//   diel.set_region(40, 60, 40, 60, 40, 60);
//   diel.set_material_id(3);

//   fdtd.add_geometry(&all);
//   fdtd.add_geometry(&diel);

//   // Excitation
//   Gaussm gm;
//   gm.set_parameters(1, 500e12, 300e12);

//   BartlettExcitation ex(&gm);
//   ex.set_soft(false);
//   ex.set_region(50, 50, 50, 50, 50, 50);
//   ex.set_polarization(0.0, 1.0, 0.0);
  
//   fdtd.add_excitation("modgauss", &ex);

//   // Results
//   grid_point p;
//   p.x = 50;
//   p.y = 25;
//   p.z = 50;
//   PointResult res1(p);
//   PointDFTResult pdft(100e12, 600e12, 50);
//   pdft.set_point(p);

//   fdtd.add_result("res1", &res1);
//   fdtd.add_result("pdft", &pdft);

//   AsciiDataWriter adw1;
//   adw1.set_filename("t_field_50.txt");

//   AsciiDataWriter adw2;
//   adw2.set_filename("t_field_dft_50.txt");

//   fdtd.add_datawriter("adw1", &adw1);
//   fdtd.add_datawriter("adw2", &adw2);
//   fdtd.map_result_to_datawriter("res1", "adw1");
//   fdtd.map_result_to_datawriter("pdft", "adw2");

//   NetCDFDataWriter ncdw;
//   ncdw.set_filename("yz_plane.nc");

//   fdtd.add_datawriter("ncdw", &ncdw);

//   PlaneResult pr1;
//   pr1.set_name("yzplane");
//   pr1.set_plane(p, BACK);

//   fdtd.add_result("pr1", &pr1);
//   fdtd.map_result_to_datawriter("pr1", "ncdw");

//   SourceDFTResult sdftr(gm, 100e12, 600e12, 50);
//   sdftr.set_time_param(0, 100, 0);
//   fdtd.add_result("sdftr", &sdftr);

//   AsciiDataWriter adw5;
//   adw5.set_filename("src_dft.txt");
//   fdtd.add_datawriter("adw5", &adw5);
  
//   fdtd.map_result_to_datawriter("sdftr", "adw5");

//   fdtd.set_time_steps(100);
//   fdtd.run(rank, size);
// }

// void pml_test(int rank, int size)
// {
//   FDTD fdtd;
  
//   fdtd.set_grid_size(900, 317, 79);

//   fdtd.set_grid_deltas(18.75e-9, 18.75e-9, 18.75e-9);
//   //fdtd.set_time_delta(3.1250e-17);

//    Pml front(VP, 1.0), back(VP, 1.0), left(VP, 1.0), right(VP, 1.0),
//      top(VP, 1.0), bottom(VP, 1.0);
//    front.set_thickness(4);
//    back.set_thickness(4);
//    left.set_thickness(4);
//    right.set_thickness(4);
//    top.set_thickness(4);
//    bottom.set_thickness(4);

//   Ewall ewall;
//   //UPml front, back, left, right, top, bottom;
//   //front.set_thickness(4); back.set_thickness(4);
//   //left.set_thickness(4); right.set_thickness(4);
//   //top.set_thickness(4); bottom.set_thickness(4);

//   fdtd.set_boundary(FRONT, &ewall);
//   fdtd.set_boundary(BACK, &ewall);
//   fdtd.set_boundary(BOTTOM, &ewall);
//   fdtd.set_boundary(TOP, &ewall);
//   fdtd.set_boundary(LEFT, &ewall);
//   fdtd.set_boundary(RIGHT, &ewall);

// //    fdtd.set_boundary(FRONT, &front);
// //    fdtd.set_boundary(BACK, &back);
// //    fdtd.set_boundary(BOTTOM, &bottom);
// //    fdtd.set_boundary(TOP, &top);
// //    fdtd.set_boundary(LEFT, &left);
// //    fdtd.set_boundary(RIGHT, &right);


//   MaterialLib mats; 
//   Material mat; // defaults to free space
//   mats.add_material(mat);
  
//   // The library stores copies. 
//   mat.set_epsilon(2.2);
//   mat.set_name("Substrate");
//   mats.add_material(mat);

//   mat.set_epsilon(1);
//   mat.set_name("Silver");
//   //mat.set_collision_freq(57e12); // THz
//   //mat.set_plasma_freq(2 * PI * 2000e+12); // THz * 2 * pi
//   mat.set_collision_freq(1.4e+14);
//   mat.set_plasma_freq(2 * PI * 1.85e+15);
//   //mats.add_material(mat);

//   fdtd.load_materials(mats);

//   // Global coordinates. 
//   Box all;
//   all.set_region(0, 100, 0, 50, 0, 50);
//   all.set_material_id(1);

//   Box metal1;
//   metal1.set_region(40, 65, 10, 40, 10, 40); // UNSTABLE
//   metal1.set_material_id(3);

// //    Box metal2;
// //    metal2.set_region(45, 50, 5, 46, 35, 60);
// //    metal2.set_material_id(3);

// //    Sphere sp1;
// //    sp1.set_centre(grid_point(40, 10, 10));
// //    sp1.set_radius(5);
// //    sp1.set_material_id(3);

//   fdtd.add_geometry(&all);
//   //fdtd.add_geometry(&metal1);
//   //fdtd.add_geometry(&sp1);
//   // fdtd.add_geometry(&metal2);

//   // Excitation
//   Gaussm gm;
//   gm.set_parameters(10, 200e12, 100e12);

//   Excitation ex(&gm);
//   //BartlettExcitation ex(gm);
//   ex.set_soft(false);
//   ex.set_region(20, 20, 25, 25, 25, 25);
//   //ex.set_region(20, 20, 6, 13, 6, 13);
//   ex.set_polarization(0.0, 1.0, 0.0);

//   fdtd.add_excitation("modgauss", &ex);

//   // Results
// //    grid_point p;
// //    p.x = 35;
// //    p.y = 25;
// //    p.z = 25;
// //    PointResult res1(p);
// //    PointDFTResult pdft(5e12, 600e12, 120);
// //    pdft.set_point(p);

// //    fdtd.add_result("res1", &res1);
// //    fdtd.add_result("pdft", &pdft);
  
// //    AsciiDataWriter adw1(rank, size);
// //    adw1.set_filename("t_field_20.txt");

// //    AsciiDataWriter adw2(rank, size);
// //    adw2.set_filename("t_field_dft_20.txt");

// //    fdtd.add_datawriter("adw1", &adw1);
// //    fdtd.add_datawriter("adw2", &adw2);
// //    fdtd.map_result_to_datawriter("res1", "adw1");
// //    fdtd.map_result_to_datawriter("pdft", "adw2");

// //    grid_point p3(35, 10, 9);
// //    PointResult pres60(p3);
// //    PointDFTResult pdft60(5e12, 600e12, 120);  

// //    fdtd.add_result("pres60", &pres60);
// //    fdtd.add_result("pdft60", &pdft60);

// //    AsciiDataWriter adwp60(rank, size);
// //    adwp60.set_filename("t_field_21.txt");

// //    AsciiDataWriter adwp60dft(rank, size);
// //    adwp60dft.set_filename("t_field_dft_21.txt");

// //    fdtd.add_datawriter("adwp60", &adwp60);
// //    fdtd.add_datawriter("adwp60dft", &adwp60dft);
// //    fdtd.map_result_to_datawriter("pres60", "adwp60");
// //    fdtd.map_result_to_datawriter("pdft60", "adwp60dft");

// //    grid_point p2;
// //    p2.x = 4;
// //    p2.y = 10;
// //    p2.z = 10;
// //    PointResult res4(p2);
// //    PointDFTResult p2dft(100e12, 600e12, 50);
// //    p2dft.set_point(p2);

// //    fdtd.add_result("res4", &res4);
// //    fdtd.add_result("p2dft", &p2dft);
  
// //    AsciiDataWriter adw12(rank, size);
// //    adw12.set_filename("t_field_4.txt");

// //    AsciiDataWriter adw13(rank, size);
// //    adw13.set_filename("t_field_dft_4.txt");

// //    fdtd.add_datawriter("adw12", &adw12);
// //    fdtd.add_datawriter("adw13", &adw13);
// //    fdtd.map_result_to_datawriter("res4", "adw12");
// //    fdtd.map_result_to_datawriter("p2dft", "adw13");

// //    MatlabDataWriter mdw(rank, size);
// //    mdw.set_filename("test.mat");
// //    fdtd.add_datawriter("mdw", &mdw);
// //    //fdtd.map_result_to_datawriter("pres60", "mdw");

// //    FarfieldResult farfield; 
// //    farfield.set_mpi_rank_size(rank, size);
// //    farfield.set_angles(0, 90, 90, 90, 45);
// //    farfield.set_freq_start(300e12);
// //    farfield.set_freq_stop(300e12);
// //    farfield.set_num_freq(1);
// //    region_t h;
// //    h.xmin = 20; h.xmax = 60; h.ymin = 5; h.ymax = 15; 
// //    h.zmin = 5; h.zmax = 15; 
// //    farfield.set_huygens(h);
// //    farfield.set_time_param(0, 49, 0);

//   //fdtd.add_result("farfield", &farfield);
//   //fdtd.map_result_to_datawriter("farfield", "mdw");

// //    NetCDFDataWriter ncdw(rank, size);
// //    ncdw.set_filename("yz_plane.nc");

// //    fdtd.add_datawriter("ncdw", &ncdw);

// //    PlaneResult pr1;
// //    pr1.set_name("ex-yzplane4");
// //    pr1.set_plane(grid_point(4, 10, 10), FRONT);
// //    pr1.set_field(FC_EX);
  
// //    PlaneResult pr2;
// //    pr2.set_name("ey-xzplane");
// //    pr2.set_plane(p, LEFT);
// //    pr2.set_field(FC_EY);

// //     PlaneResult pr3;
// //     pr3.set_name("ey-yzplane4");
// //     pr3.set_plane(grid_point(4, 10, 10), FRONT);
// //     pr3.set_field(FC_EY);

// //     PlaneResult pr4;
// //     pr4.set_name("ey-yzplane15");
// //     pr4.set_plane(grid_point(15, 10, 10), FRONT);
// //     pr4.set_field(FC_EZ);

// //     PlaneResult pr5;
// //     pr5.set_name("ey-xyplane5");
// //     pr5.set_plane(p, BOTTOM);
// //     pr5.set_field(FC_EY);

// //     fdtd.add_result("pr1", &pr1);
// //     fdtd.add_result("pr2", &pr2);
// //     fdtd.add_result("pr3", &pr3);
// //     fdtd.add_result("pr4", &pr4);
// //     fdtd.add_result("pr5", &pr5);

// //     fdtd.map_result_to_datawriter("pr1", "ncdw");
// //     fdtd.map_result_to_datawriter("pr2", "ncdw");
// //     fdtd.map_result_to_datawriter("pr3", "ncdw");
// //     fdtd.map_result_to_datawriter("pr4", "ncdw");
// //     fdtd.map_result_to_datawriter("pr5", "ncdw");

// //     SourceDFTResult sdftr(gm, 100e12, 600e12, 50);
// //     sdftr.set_time_param(0, 500, 0);
// //     fdtd.add_result("sdftr", &sdftr);

// //     AsciiDataWriter adw5(rank, size);
// //     adw5.set_filename("src_dft.txt");
// //     fdtd.add_datawriter("adw5", &adw5);

// //     fdtd.map_result_to_datawriter("sdftr", "adw5");

// //     SourceTimeResult srctr(gm);
// //     AsciiDataWriter adw8(rank, size);
// //     adw8.set_filename("src.txt");

// //     fdtd.add_result("srctr", &srctr);
// //     fdtd.add_datawriter("adw8", &adw8);

// //     fdtd.map_result_to_datawriter("srctr", "adw8");

//    fdtd.set_time_steps(100);

   
// #ifdef USE_OPENMP
//    // Test the OpenMP
//    time_t start, now;
//    clock_t cpu_start, cpu_now;

//    for (int numthreads = 1; numthreads <= omp_get_max_threads(); numthreads++)
//      {
//        omp_set_num_threads(numthreads);
//        start = time(NULL);
//        cpu_start = clock();
//        fdtd.run(rank, size);
//        now = time(NULL);
//        cpu_now = clock();
       
//        cout << numthreads << " of " 
// 	    << omp_get_max_threads() << " threads took " 
// 	    << now - start << " wall clock seconds, and "
// 	    << (cpu_now - cpu_start) / static_cast<double>(CLOCKS_PER_SEC)
// 	    << " cpu seconds." << endl;
//      }

// #else
//    fdtd.run(rank, size);
// #endif
// }

// void takakura_test(int rank, int size)
// {
//   FDTD fdtd;
  
//   fdtd.set_grid_size(50, 320, 266);
//   fdtd.set_grid_deltas(18.75e-9, 18.75e-9, 18.75e-9);
//   //fdtd.set_time_delta(3.1250e-17);

//   Pml front(VP, 1.0), back(VP, 1.0), left(VP, 1.0), right(VP, 1.0),
//     top(VP, 1.0), bottom(VP, 1.0);
//   front.set_thickness(5);
//   back.set_thickness(5);
//   left.set_thickness(5);
//   right.set_thickness(5);
//   top.set_thickness(5);
//   bottom.set_thickness(5);

//   fdtd.set_boundary(FRONT, &front);
//   fdtd.set_boundary(BACK, &back);
//   fdtd.set_boundary(BOTTOM, &bottom);
//   fdtd.set_boundary(TOP, &top);
//   fdtd.set_boundary(LEFT, &left);
//   fdtd.set_boundary(RIGHT, &right);

//   MaterialLib mats; 
//   Material mat; // defaults to free space
//   mats.add_material(mat);
  
//   // The library stores copies. 
//   mat.set_epsilon(2.2);
//   mat.set_name("Substrate");
//   mats.add_material(mat);

//   mat.set_epsilon(1);
//   mat.set_name("Silver");
//   //mat.set_collision_freq(57e12); // THz
//   //mat.set_plasma_freq(2 * PI * 2000e+12); // THz * 2 * pi
//   mat.set_collision_freq(1.4e+14);
//   mat.set_plasma_freq(2 * PI * 1.85e+15);
//   mats.add_material(mat);

//   fdtd.load_materials(mats);

//   // Global coordinates. 
//   Box all;
//   all.set_region(0, 50, 0, 320, 0, 266);
//   all.set_material_id(1);

//   Box metal1;
//   metal1.set_region(6, 43, 6, 155, 106, 160); // UNSTABLE
//   metal1.set_material_id(3);

//   Box metal2;
//   metal2.set_region(6, 43, 164, 313, 106, 160);
//   metal2.set_material_id(3);

//   fdtd.add_geometry(&all);
//   fdtd.add_geometry(&metal1);
//   fdtd.add_geometry(&metal2);

//   // Excitation
//   Gaussm gm;
//   gm.set_parameters(10, 500e12, 300e12);

//   //Excitation ex(&gm);
//   BartlettExcitation ex(&gm);
//   ex.set_soft(true);
//   ex.set_region(6, 43, 6, 313, 15, 15);
//   ex.set_polarization(0.0, 1.0, 0.0);

//   fdtd.add_excitation("modgauss", &ex);

//   // Results
//   grid_point p;
//   p.x = 24;
//   p.y = 159;
//   p.z = 75;
//   PointResult res1(p);
//   PointDFTResult pdft(5e12, 600e12, 120);
//   pdft.set_point(p);

//   fdtd.add_result("res1", &res1);
//   fdtd.add_result("pdft", &pdft);
  
//   AsciiDataWriter adw1;
//   adw1.set_filename("t_field_75.txt");

//   AsciiDataWriter adw2;
//   adw2.set_filename("t_field_dft_75.txt");

//   fdtd.add_datawriter("adw1", &adw1);
//   fdtd.add_datawriter("adw2", &adw2);
//   fdtd.map_result_to_datawriter("res1", "adw1");
//   fdtd.map_result_to_datawriter("pdft", "adw2");

//   grid_point p3(24, 159, 200);
//   PointResult pres60(p3);
//   PointDFTResult pdft60(5e12, 600e12, 120);  

//   fdtd.add_result("pres60", &pres60);
//   fdtd.add_result("pdft60", &pdft60);

//   AsciiDataWriter adwp60;
//   adwp60.set_filename("t_field_200.txt");

//   AsciiDataWriter adwp60dft;
//   adwp60dft.set_filename("t_field_dft_200.txt");

//   fdtd.add_datawriter("adwp60", &adwp60);
//   fdtd.add_datawriter("adwp60dft", &adwp60dft);
//   fdtd.map_result_to_datawriter("pres60", "adwp60");
//   fdtd.map_result_to_datawriter("pdft60", "adwp60dft");

//   NetCDFDataWriter ncdw;
//   ncdw.set_filename("takakura-test.nc");

//   fdtd.add_datawriter("ncdw", &ncdw);

//   PlaneResult pr1;
//   pr1.set_name("ex-yzplane4");
//   pr1.set_plane(grid_point(24, 159, 75), FRONT);
//   pr1.set_field(FC_EY);
//   pr1.set_time_param(1, 5000, 25);
 
//   fdtd.add_result("pr1", &pr1);
//   fdtd.map_result_to_datawriter("pr1", "ncdw");

//   fdtd.set_time_steps(2000);
//   fdtd.run(rank, size);
// }

// void laser_test(int rank, int size)
// {
//   FDTD fdtd;
  
//   fdtd.set_grid_size(75, 21, 21);
//   fdtd.set_grid_deltas(18.75e-9, 18.75e-9, 18.75e-9);
//   //fdtd.set_time_delta(3.1250e-17);

//   Pml front(VP, 1.0), back(VP, 1.0), left(VP, 1.0), right(VP, 1.0),
//     top(VP, 1.0), bottom(VP, 1.0);
//   front.set_thickness(4);
//   back.set_thickness(4);
//   left.set_thickness(4);
//   right.set_thickness(4);
//   top.set_thickness(4);
//   bottom.set_thickness(4);

//   fdtd.set_boundary(FRONT, &front);
//   fdtd.set_boundary(BACK, &back);
//   fdtd.set_boundary(BOTTOM, &bottom);
//   fdtd.set_boundary(TOP, &top);
//   fdtd.set_boundary(LEFT, &left);
//   fdtd.set_boundary(RIGHT, &right);

//   MaterialLib mats; 
//   Material mat; // defaults to free space
//   mats.add_material(mat);
  
//   // The library stores copies. 
//   mat.set_epsilon(2.2);
//   mat.set_name("Substrate");
//   mats.add_material(mat);

//   mat.set_epsilon(1);
//   mat.set_name("Silver");
//   //mat.set_collision_freq(57e12); // THz
//   //mat.set_plasma_freq(2 * PI * 2000e+12); // THz * 2 * pi
//   mat.set_collision_freq(1.4e+14);
//   mat.set_plasma_freq(2 * PI * 1.85e+15);
//   mats.add_material(mat);

//   fdtd.load_materials(mats);

//   // Global coordinates. 
//   Box all;
//   all.set_region(0, 75, 0, 21, 0, 21);
//   all.set_material_id(1);

//   Box metal1;
//   metal1.set_region(40, 65, 5, 14, 5, 14); // UNSTABLE
//   metal1.set_material_id(3);
// }



// void coupler_test(int rank, int size)
// {
//   FDTD fdtd;
  
//   field_t start_f = 12e9;
//   field_t stop_f = 18e9;
//   field_t centre_f = 15e9;
//   field_t d_f = 3e9;

//   unsigned int xlen = 600, ylen = 158 * 2 + 1, zlen = 79;
//   float dx = 0.05e-3, dy = 0.1e-3, dz = 0.1e-3;

//   fdtd.set_grid_size(xlen, ylen, zlen);
//   fdtd.set_grid_deltas(dx, dy, dz);

//   //fdtd.set_time_delta(3.1250e-17);

//   cout << "Ku Band waveguide aperture coupler test. Grid is "
//        << xlen << "x" << ylen << "x" << zlen << ", "
//        << "\ntime step size is " << fdtd.get_time_delta() << "." << endl;

//   Pml front(VP, 1.0), back(VP, 1.0);
//   front.set_thickness(4);
//   back.set_thickness(4);

//   Ewall ewall;

//   fdtd.set_boundary(FRONT, &front);
//   fdtd.set_boundary(BACK, &back);
//   fdtd.set_boundary(BOTTOM, &ewall);
//   fdtd.set_boundary(TOP, &ewall);
//   fdtd.set_boundary(LEFT, &ewall);
//   fdtd.set_boundary(RIGHT, &ewall);

//   MaterialLib mats; 
//   Material mat; // defaults to free space
//   mats.add_material(mat);
  
//   fdtd.load_materials(mats);

//   // Global coordinates. 
//   Box all;
//   all.set_region(0, xlen, 0, ylen, 0, zlen);
//   all.set_material_id(1);

//   Box metal1;
//   metal1.set_region(0, xlen, 159, 160, 0, zlen); 
//   metal1.set_material_id(0);

//   fdtd.add_geometry(&all);
//   fdtd.add_geometry(&metal1);

//   unsigned int A = 0, B = 0;

//   unsigned int L[] = {static_cast<unsigned int>(3e-3 / dx), 
//                       static_cast<unsigned int>(2.24e-3 / dx),
//                       static_cast<unsigned int>(2.19e-3 / dx),
//                       static_cast<unsigned int>(2.43e-3 / dx),
//                       static_cast<unsigned int>(2.19e-3 / dx),
//                       static_cast<unsigned int>(2.24e-3 / dx)};

//   unsigned int S[] = {static_cast<unsigned int>(1.52e-3 / dx),
//                       static_cast<unsigned int>(2.01e-3 / dx),
//                       static_cast<unsigned int>(3.07e-3 / dx),
//                       static_cast<unsigned int>(3.07e-3 / dx),
//                       static_cast<unsigned int>(2.01e-3 / dx),
//                       static_cast<unsigned int>(1.52e-3 / dx)};

//   Box apertures[6];

//   for (int idx = 0; idx < 6; idx++)
//   {
//     A += L[idx];
//     B = A + S[idx];

//     apertures[idx].set_region(A, B, 159, 160, 0, zlen);
//     apertures[idx].set_material_id(1);

//     cout << "Adding aperture, x = " << A << " to " << B << endl;

//     A += S[idx];

//     fdtd.add_geometry(&apertures[idx]);
//   }

//   // Excitation
//   Gaussm gm;
//   gm.set_parameters(10, d_f, centre_f);

//   WaveguideExcitation ex(&gm);
//   ex.set_soft(true);
//   ex.set_region(10, 11, 0, ylen, 0, zlen);
//   ex.set_polarization(0.0, 0.0, 1.0);
//   ex.set_mode(0, 1, 0);

//   fdtd.add_excitation("modgauss", &ex);

//   // Results
//   grid_point p;
//   p.x = xlen / 2;
//   p.y = ylen / 2;
//   p.z = zlen / 2;
  
//   MatlabDataWriter mdw;
//   mdw.set_filename("coupler.mat");
//   fdtd.add_datawriter("mdw", &mdw);

// //  NetCDFDataWriter ncdw(rank, size);
// //  ncdw.set_filename("coupler.nc");
// //  fdtd.add_datawriter("ncdw", &ncdw);

// //  PlaneResult pr1;
// //  pr1.set_name("ez_xyplane");
// //  pr1.set_plane(p, TOP);
// //  pr1.set_field(FC_EZ);
 
// //  fdtd.add_result("ez_xyplane", &pr1);
// //  fdtd.map_result_to_datawriter("ez_xyplane", "ncdw");

//   // S Parameters
//   PowerResult s11(start_f, stop_f, 120), s12(start_f, stop_f, 120),
//     s13(start_f, stop_f, 120), s14(start_f, stop_f, 120);

//   region_t s11r, s12r, s13r, s14r;

//   s11r.xmin = 20; s11r.xmax = 20; s11r.ymin = 0; s11r.ymax = 158;
//   s11r.zmin = 0; s11r.zmax = zlen;

//   s12r = s13r = s14r = s11r; 

//   s12r.xmin = xlen - 20; s12r.xmax = xlen - 20;
//   s13r.xmin = xlen - 20; s13r.xmax = xlen - 20;
//   s13r.ymin = 160; s13r.ymax = ylen;
//   s14r.ymin = 160; s14r.ymax = ylen;

//   s11.set_region(s11r);
//   s12.set_region(s12r);
//   s13.set_region(s13r);
//   s14.set_region(s14r);

//   fdtd.add_result("s11", &s11);
//   fdtd.add_result("s12", &s12);
//   fdtd.add_result("s13", &s13);
//   fdtd.add_result("s14", &s14);

//   fdtd.map_result_to_datawriter("s11", "mdw");
//   fdtd.map_result_to_datawriter("s12", "mdw");
//   fdtd.map_result_to_datawriter("s13", "mdw");
//   fdtd.map_result_to_datawriter("s14", "mdw");

//   // The source, just for fun. 
//   SourceDFTResult sdftr(gm, 100e12, 600e12, 50);
//   fdtd.add_result("sdftr", &sdftr);
//   fdtd.map_result_to_datawriter("sdftr", "mdw");

//   SourceTimeResult srctr(gm);
//   fdtd.add_result("src", &srctr);
//   fdtd.map_result_to_datawriter("src", "mdw");


//   // Let's begin.
//   fdtd.set_time_steps(10000);
//   fdtd.run(rank, size);
// }
