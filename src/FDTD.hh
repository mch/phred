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

#ifndef FDTD_H
#define FDTD_H

#include <map>
#include <vector>
#include <string>

using namespace std;

#include <boost/shared_ptr.hpp>

using namespace boost;

#include "Types.hh"
#include "Boundaries/BoundaryCondition.hh"
#include "MaterialLib.hh"
#include "Excitations/Excitation.hh"
#include "Results/Result.hh"
#include "DataWriters/DataWriter.hh"
#include "Constants.hh"
#include "Exceptions.hh"
#include "SimpleSDAlg.hh"
#include "ProblemGeometry.hh"
#include "FreqGrid.hh"

/**
 * This is a sort of convience wrapper object that runs the
 * simulation. It is primarily intended to be used from python. It
 * holds lists of results, excitations, data writers, etc, manages
 * thier memory and life cycles, etc. 
 */
class FDTD
{
public:
  FDTD();
  virtual ~FDTD();

  /**
   * Set the grid material
   */
  void set_grid_material(const char *material);

  /**
   * Set the number of time steps to run the simulation for. 
   */ 
  void set_time_steps(unsigned int t);

  /**
   * Set the length of time the simulation should run for. 
   * The number of time steps is automatically calculated.
   */ 
  void set_time(float t);

  /**
   * Set the global (entire problem) grid size, in meters. 
   */
  void set_grid_size(float x, float y, float z);

  /**
   * Set the center of the grid. Defaults to (0,0,0).
   */ 
  void set_grid_centre(float x, float y, float z);

  /**
   * Returns the grid center as a point.
   */ 
  inline point get_grid_centre()
  { return geometry_.get_grid_centre(); }

  /**
   * Returns the grid size as a point representing the length along
   * each axis.
   */ 
  inline point get_grid_size()
  { return geometry_.get_grid_size(); }

  /**
   * Set the grid deltas, the size of the cells. The time delta is
   * automatically computed using the stability condition. 
   */
  void set_grid_deltas(field_t dx, field_t dy, field_t dz);

  /**
   * Set a time delta different from the one computed by
   * set_grid_deltas(). 
   */
  void set_time_delta(field_t dt);

  /**
   * Returns the time delta. This is computed when set_grid_deltas is
   * called.
   */
  field_t get_time_delta();

  /**
   * Returns the number of cells in the x direction
   */ 
  unsigned int get_num_x_cells();

  /**
   * Returns the number of cells in the y direction
   */ 
  unsigned int get_num_y_cells();

  /**
   * Returns the number of cells in the z direction
   */ 
  unsigned int get_num_z_cells();

  /**
   * Set a boundary condition object. A pointer to the object is
   * stored, so don't, for the love of god and all that is holy,
   * delete the object or allow it go go out of scope before this
   * object does.
   *
   * This now uses a shared_ptr, so its ok. 
   */
  void set_boundary(Face face, shared_ptr<BoundaryCond> bc);

  /**
   * Load a material library to use.
   */
  void load_materials(shared_ptr<MaterialLib> matlib);

  /**
   * Let us know about an excitation, or replace one of the
   * same name.  You are passing in a pointer; don't dispose of the
   * object!
   */
  void add_excitation(const char *name, shared_ptr<Excitation> ex);

  /**
   * Add a result object, or replace one of the same name. 
   */
  void add_result(const char *name, shared_ptr<Result> r);

  /** 
   * Add a datawriter, or replace one of the same name. 
   */
  void add_datawriter(const char *name, shared_ptr<DataWriter> dw);

  /**
   * Add a CSG object to the grid
   */
  void add_object(string material, shared_ptr<CSGObject> obj);

  /**
   * Map a results to a DataWriter. Some DataWriters cannot accept
   * more than one result, so this may throw and exception. Use this
   * function after the results and datawriters have been added. 
   *
   * @param result the name of the result
   * @param dw the name of the data writer
   */
  void map_result_to_datawriter(const char *result, const char *dw);

  /**
   * Run the simulation for N time steps. 
   */
  void run();

  /**
   * This is to allow for testing, so that indiviual components can be
   * tested from Python seperatly.
   */ 
  inline shared_ptr<Grid> get_grid()
  { return grid_; }

protected:
  /**
   * The grid to operate on; what kind of grid specifically is decided
   * at run time from material properties.
   */
  shared_ptr<Grid> grid_;

  /**
   * Global grid information
   */
  GridInfo global_ginfo_;

  /**
   * Local grid information, which actually corresponds to
   * grid_. Created by a domain decomposition algorithm. 
   */
  GridInfo local_ginfo_;

  /**
   * Our E excitations 
   */
  map<string, shared_ptr<Excitation> > e_excitations_;

  /**
   * Our H excitations 
   */
  map<string, shared_ptr<Excitation> > h_excitations_;

  /**
   * Our results
   */
  map<string, shared_ptr<Result> > results_;

  /**
   * Our data writers
   */
  map<string, shared_ptr<DataWriter> > datawriters_;

  /**
   * A map that tells which result goes to which data writer. Results
   * may go to multiple data writers, and some data writers can
   * recieve multiple results. 
   */
  vector< pair<string, string> > r_dw_map_;

  /**
   * Geometry objects
   */
  ProblemGeometry geometry_;

  /**
   * Material library
   */
  shared_ptr<MaterialLib> mlib_;

  /**
   * Number of time steps to go for. 
   */ 
  unsigned int time_steps_; 

  /** 
   * Call LifeCycle::init() 
   */
  //template<class T, class A>
  //void init_objs();

  /** 
   * Call LifeCycle::deinit() 
   */
  //template<class T, class A>
  //void deinit_objs();  

  /**
   * Adds the results to the datawriters. This is called from run()
   * after the datawriters and results have been initialized since
   * some DataWriters need to be initialized before variables can be
   * added to them. 
   */
  void setup_datawriters();

};

#endif // FDTD_H
