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
   * The grid to operate on
   */
  Grid grid_;

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
   * Our excitations
   */
  map<string, Excitation *> excitations_;

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
   * Let us know about an excitation, or replace one of the same name.
   * You are passing in a pointer; don't dispose of the object!
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
   * @param steps number of timesteps to run for. 
   */
  void run(int rank, int size, unsigned int steps);

};

#endif // FDTD_H
