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

#ifndef FDTD_H
#define FDTD_H

#include <map>
#include <vector>
#include <string>

using namespace std;

#include "Types.hh"
#include "BoundaryCondition.hh"
#include "MaterialLib.hh"
#include "Excitation.hh"
#include "Result.hh"
#include "DataWriter.hh"
#include "Constants.hh"
#include "Exceptions.hh"
#include "SimpleSDAlg.hh"
#include "Geometry.hh"
#include "FreqGrid.hh"

/**
 * This is a sort of convience wrapper object that runs the
 * simulation. It is primarily intended to be used from python. It
 * holds lists of results, excitations, data writers, etc, manages
 * thier memory and life cycles, etc. 
 *
 */
class FDTD
{
private:
protected:
  /**
   * The grid to operate on; what kind of grid specifically is decided
   * at run time from material properties.
   */
  Grid *grid_;

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
  map<string, Excitation *> e_excitations_;

  /**
   * Our H excitations 
   */
  map<string, Excitation *> h_excitations_;

  /**
   * Our results
   */
  map<string, Result *> results_;

  /**
   * Our data writers
   */
  map<string, DataWriter *> datawriters_;

  /**
   * A map that tells which result goes to which data writer. Results
   * may go to multiple data writers, and some data writers can
   * recieve multiple results. 
   */
  vector< pair<string, string> > r_dw_map_;

  /**
   * Geometry objects
   */
  vector<Geometry *> geometry_;

  /**
   * Material library
   */
  MaterialLib *mlib_;

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

public:
  FDTD();
  virtual ~FDTD();

  /**
   * Set the number of time steps to run the simulation for. 
   */ 
  void set_time_steps(unsigned int t);

  /**
   * Set the global (entire problem) grid size
   */
  void set_grid_size(unsigned int x, 
                     unsigned int y, unsigned int z);

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
   * Set a boundary condition object. A pointer to the object is
   * stored, so don't, for the love of god and all that is holy,
   * delete the object or allow it go go out of scope before this
   * object does.
   */
  void set_boundary(Face face, BoundaryCond *bc);

  /**
   * Load a material library to use.
   */
  void load_materials(MaterialLib &matlib);

  /**
   * Let us know about an excitation, or replace one of the
   * same name.  You are passing in a pointer; don't dispose of the
   * object!
   */
  void add_excitation(const char *name, Excitation *ex);

  /**
   * Add a result object, or replace one of the same name. 
   */
  void add_result(const char *name, Result *r);

  /** 
   * Add a datawriter, or replace one of the same name. 
   */
  void add_datawriter(const char *name, DataWriter *dw);

  /**
   * Add a geometry object to the grid
   */
  void add_geometry(Geometry *geom);

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
   *
   * @param rank the process rank in MPI
   * @param size the number of ranks in the MPI communicator
   */
  void run(int rank, int size);

  /**
   * Looks through the geometry stack and returns a reference to the
   * geometry that a particular point belongs to. 
   *
   * @param p the point to find the geometry for
   */
  Geometry &find_geometry(unsigned int x,
                          unsigned int y,
                          unsigned int z);

};

#endif // FDTD_H
