/* 
   Phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004-2005 Matt Hughes <mhughe@uvic.ca>

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

#include "MetaFDTD.hh"
#include "Globals.hh"

// The Meta part:
#include "GridUpdate.hh"


MetaFDTD::MetaFDTD()
  : FDTD(), mt_(METAFDTD_THREE)
{}

MetaFDTD::~MetaFDTD()
{}

// This is a straight copy from FDTD.cc. This really needs some
// refactoring.

// CHOP THIS UP; make helper functions
void MetaFDTD::run()
{
  // Check that we have every thing we need...
  // 1) Material Library
  // 2) ...

  // Determine the number of cells in the Grid from the size of the
  // grid box and the maximum excitation frequency, unless overridden.
  //...

  // Subdivide the grid using a domain decomosition algorithm.
  if (!quiet)
    cout << "Performing domain decomposition..." << endl;

  SubdomainAlg *alg = 0;
  alg = new MPISubdomainAlg();
  
  local_ginfo_ = alg->decompose_domain(global_ginfo_);

  delete alg;

  if (!quiet)
  {
    cout
      << "Size of the entire (global) grid: \n\t-> " 
      << local_ginfo_.global_dimx_ << " x " 
      << local_ginfo_.global_dimy_ << " x " 
      << local_ginfo_.global_dimz_ << ".\n"
    
      << "Size of the local grid on rank " << MPI_RANK << ":\n\t-> "
      << local_ginfo_.dimx_ << " x " 
      << local_ginfo_.dimy_ << " x " 
      << local_ginfo_.dimz_ << ".\n"

      << "Size of the local grid excluding ghost cells:\n\t-> "
      << local_ginfo_.dimx_no_sd_ << " x " 
      << local_ginfo_.dimy_no_sd_ << " x " 
      << local_ginfo_.dimz_no_sd_ << ".\n"

      << "Within the global grid, the local grid on rank " 
      << MPI_RANK << " starts at:\n\t-> "
      << local_ginfo_.start_x_ << " x " 
      << local_ginfo_.start_y_ << " x " 
      << local_ginfo_.start_z_ << ".\n"

    << "Within the global grid, the local grid on rank " 
      << MPI_RANK << ", excluding ghost cells, \nstarts at:\n\t-> "
      << local_ginfo_.start_x_no_sd_ << " x " 
      << local_ginfo_.start_y_no_sd_ << " x " 
      << local_ginfo_.start_z_no_sd_ << ".\n\n";
  }

  // Decide what grid to used from materials
  bool freqgrid = false;
  map<string, Material>::const_iterator miter = mlib_->get_material_iter_begin();
  map<string, Material>::const_iterator miter_e = mlib_->get_material_iter_end();
  
  while(miter != miter_e)
  {
    if (((*miter).second).get_collision_freq() != 0.0 ||
        ((*miter).second).get_plasma_freq() != 0.0)
    {
      freqgrid = true; 
      break;
    }
    miter++;
  }

  if (freqgrid)
    {
      grid_ = shared_ptr<Grid>(new FreqGrid());
#ifdef DEBUG
      cout << "Using freq grid. " << endl;
#endif
    }
  else
    {
#ifdef DEBUG
      cout << "Using simple grid. " << endl;
#endif
      grid_ = shared_ptr<Grid>(new Grid());
    }

  if (!quiet)
    cout << "Initializing grid..." << endl;

  grid_->setup_grid(local_ginfo_);

  grid_->load_materials(mlib_);

  geometry_.init(*grid_);
  grid_->load_geometry(&geometry_);

  grid_->set_define_mode(false);

  if (!quiet)
  {
    cout << "Grid is " << get_num_x_cells() << "x" << get_num_y_cells()
         << "x" << get_num_z_cells() << " cells in size.\n";

    if (MPI_SIZE > 1)
    {
      cout << "The Grid on this node is " 
           << local_ginfo_.dimx_ << "x" 
           << local_ginfo_.dimy_ << "x" 
           << local_ginfo_.dimz_ << " cells in size.\n";
    }

    cout << "Time step size is " << get_time_delta() 
         << ". Executing for " << time_steps_ << " time steps.\n\n";
    cout << "Initializing results, data writers, and excitations..." << endl;
  }
  // Life cycle init
  map<string, shared_ptr<Excitation> >::iterator 
    e_eiter_b = e_excitations_.begin();
  map<string, shared_ptr<Excitation> >::iterator e_eiter = e_eiter_b;
  map<string, shared_ptr<Excitation> >::iterator 
    e_eiter_e = e_excitations_.end();
  
  if (!quiet)
    cout << "Electric excitations: \n";

  while (e_eiter != e_eiter_e) 
  {
    (*e_eiter).second->init(*grid_);
    
    if (!quiet)
      cout << "\t -> " << *((*e_eiter).second) << "\n";

    ++e_eiter;
  }

  map<string, shared_ptr<Excitation> >::iterator 
    h_eiter_b = h_excitations_.begin();
  map<string, shared_ptr<Excitation> >::iterator h_eiter = h_eiter_b;
  map<string, shared_ptr<Excitation> >::iterator 
    h_eiter_e = h_excitations_.end();
  
  if (!quiet)
    cout << "\nMagnetic excitations: \n";

  while (h_eiter != h_eiter_e) 
  {
    (*h_eiter).second->init(*grid_);

    if (!quiet)
      cout << "\t -> " << *((*h_eiter).second) << "\n";

    ++h_eiter;
  }

  if (e_eiter_b == e_eiter_e && h_eiter_b == h_eiter_e)
  {
    cout << "WARNING: NO EXCITATIONS DEFINED!\n";
  }


  map<string, shared_ptr<Result> >::iterator riter_b = results_.begin();
  map<string, shared_ptr<Result> >::iterator riter = riter_b;
  map<string, shared_ptr<Result> >::iterator riter_e = results_.end();
  
  if (!quiet)
    cout << "\nResult generators: \n";

  if (riter == riter_e)
    cout << "WARNING: NO RESULTS DEFINED!\n";

  while (riter != riter_e) 
  {
    if ((*riter).second->get_time_stop() == ~0)
      (*riter).second->set_time_stop(time_steps_);

    (*riter).second->init(*grid_);

    if (!quiet)
      cout << "\t -> " << (*(*riter).second) << "\n";

    ++riter;
  }

  map<string, shared_ptr<DataWriter> >::iterator
    dwiter_b = datawriters_.begin();
  map<string, shared_ptr<DataWriter> >::iterator dwiter = dwiter_b;
  map<string, shared_ptr<DataWriter> >::iterator
    dwiter_e = datawriters_.end();
  
  if (!quiet)
    cout << "\nDataWriters: \n";

  while (dwiter != dwiter_e) 
  {
    (*dwiter).second->init(*grid_);

    if (!quiet)
      cout << "\t -> " << (*(*dwiter).second) << "\n";

    ++dwiter;
  }

  setup_datawriters();

  // This barrier prevent ranks > 0 from continuing on if rank 0 has
  // trouble with a data writer or something. 
  MPI_Barrier(MPI_COMM_PHRED);

  if (!setup_only)
  {
    // Run
    unsigned int ts = 0;
    
    // Do data output for results that do that first thing, like GridResult
    vector< pair<string, string> >::iterator iter = r_dw_map_.begin();
    vector< pair<string, string> >::iterator iter_e = r_dw_map_.end();
    
    while (iter != iter_e)
    {
      riter = results_.find((*iter).first);
      dwiter = datawriters_.find((*iter).second);      
      
      if (riter != riter_e && dwiter != dwiter_e)
        (*dwiter).second->handle_data(0, 
                                      (*riter).second->get_pre_result(*grid_));
      ++iter;
    }
    
    // Set up the GridUpdate stuff
    GridUpdateData gud(*grid_);
    
    compute_update_regions();
    
    if (MPI_RANK == 0)
      cout << "\nStarting FDTD time stepping, running for " 
           << time_steps_ << " time steps..." << endl;

    // For optionally tracking millions of nodes per second. 
    time_t start = time(NULL);
    time_t now;

    clock_t start_cpu = clock();
    clock_t now_cpu;

    double time_total = 0;
    double time_total_cpu = 0;

    // For estimating run time
    time_t rt_start = time(NULL);
    time_t rt_now;
    unsigned int rt_steps = 9;
  
    for (ts = 1; ts <= time_steps_; ts++) {
    
      if (!extra_quiet_g && MPI_RANK == 0 && ts % 100 == 0)
      {
        rt_now = time(NULL);
        int secs = static_cast<int>((static_cast<double>(rt_now - rt_start) 
                                     / rt_steps) * (time_steps_ - ts));

        cout << "Estimated time remaining at time step " << ts << ": \n\t";
        print_elapsed_time(secs);

        rt_steps = 100;
        rt_start = time(NULL);
      }

      update_h(gud);

      if (sigterm_g) {
        break;
      }

      // Excitations
      h_eiter = h_eiter_b;
      while (h_eiter != h_eiter_e)
      {
        (*h_eiter).second->excite(*grid_, ts, H);
        ++h_eiter;
      }

      // Boundary condition application
      grid_->apply_boundaries(H);

      update_e(gud);
    
      if (sigterm_g) {
        break;
      }

      // Excitations
      e_eiter = e_eiter_b;
      while (e_eiter != e_eiter_e)
      {
        (*e_eiter).second->excite(*grid_, ts, E);
        ++e_eiter;
      }
    
      // Boundary condition application
      grid_->apply_boundaries(E);

      // Results
      iter = r_dw_map_.begin();
      iter_e = r_dw_map_.end();

      while (iter != iter_e)
      {
        riter = results_.find((*iter).first);
        dwiter = datawriters_.find((*iter).second);      
      
        if (riter != riter_e && dwiter != dwiter_e)
          (*dwiter).second->handle_data(ts, 
                                        (*riter).second->get_result(*grid_, ts));


        ++iter;
      }
    
    } // End of main loop

    now = time(NULL);
    now_cpu = clock();
    time_total = static_cast<double>(now) - static_cast<double>(start);
    time_total_cpu = static_cast<double>(now_cpu)
      - static_cast<double>(start_cpu);

    double avg_time = time_total / static_cast<double>(time_steps_);

    double avg_cpu_time = time_total_cpu / time_steps_;
    double num_mnodes = static_cast<double>(grid_->get_num_updated_nodes()) 
      / 1.0e6;

    double cpu_mnps = num_mnodes 
      / (avg_cpu_time / static_cast<double>(CLOCKS_PER_SEC));
    double wall_mnps = num_mnodes / avg_time;
    unsigned int num_nodes = grid_->get_num_updated_nodes();

    double total_mnodes, total_cpu_mnps, total_wall_mnps;
    unsigned int total_nodes;

    MPI_Reduce(&num_mnodes, &total_mnodes, 1, MPI_DOUBLE, 
               MPI_SUM, 0, MPI_COMM_PHRED);

    MPI_Reduce(&num_nodes, &total_nodes, 1, MPI_UNSIGNED, 
               MPI_SUM, 0, MPI_COMM_PHRED);

    MPI_Reduce(&cpu_mnps, &total_cpu_mnps, 1, MPI_DOUBLE, 
               MPI_SUM, 0, MPI_COMM_PHRED);

    MPI_Reduce(&wall_mnps, &total_wall_mnps, 1, MPI_DOUBLE, 
               MPI_SUM, 0, MPI_COMM_PHRED);

    if (MPI_RANK == 0)
    {
      cout << "Number of updated nodes: " << total_nodes
           << ", millions of updated nodes: " << total_mnodes << endl;
      cout << "Average wall time: " << avg_time << ", avg cpu time: "
           << avg_cpu_time / static_cast<double>(CLOCKS_PER_SEC) << endl;
      cout << "Average millions of nodes per second, w.r.t. wall clock time: " 
           << total_wall_mnps << endl;
      cout << "Average millions of nodes per second, w.r.t. CPU time: " 
           << total_cpu_mnps << endl;
      cout << "Note: MNPS w.r.t. CPU time may be incorrect. " << endl;
    }

    // Do data output for results that only produce data at the end
    iter = r_dw_map_.begin();
    iter_e = r_dw_map_.end();

    while (iter != iter_e)
    {
      riter = results_.find((*iter).first);
      dwiter = datawriters_.find((*iter).second);      
    
      if (riter != riter_e && dwiter != dwiter_e)
        (*dwiter).second->handle_data(ts, 
                                      (*riter).second->get_post_result(*grid_));
    
    
      ++iter;
    }

  } else {
    cout << "Set up only requested. Exiting. \n";
  }

  // life cycle de init
  e_eiter = e_eiter_b;
  while (e_eiter != e_eiter_e) 
  {
    (*e_eiter).second->deinit();
    ++e_eiter;
  }

  h_eiter = h_eiter_b;
  while (h_eiter != h_eiter_e) 
  {
    (*h_eiter).second->deinit();
    ++h_eiter;
  }

  riter = riter_b;
  while (riter != riter_e) 
  {
    (*riter).second->deinit();
    ++riter;
  }

  dwiter = dwiter_b;
  while (dwiter != dwiter_e) 
  {
    (*dwiter).second->deinit();
    ++dwiter;
  }

}

template<class Update>
static inline void component_correct(region_t update_r, region_t update_start,
                                     GridUpdateData &gud, Grid &grid_)
{
  // Oh man, so much code! The ranges should be in vectors so it's
  // possible to loop over them. This is super ugly. 
  region_t temp_update;

  if (update_r.xmin > update_start.xmin)
  {
    temp_update = update_start;
    temp_update.xmax = temp_update.xmin + 1;

    GridUpdateTiling<Update, GridUpdateData, PrivateGridUpdateData>
      ::grid_loop(temp_update, gud, grid_);
      
    update_start.xmin++;
  }

  if (update_r.ymin > update_start.ymin)
  {
    temp_update = update_start;
    temp_update.ymax = temp_update.ymin + 1;

    GridUpdateTiling<Update, GridUpdateData, PrivateGridUpdateData>
      ::grid_loop(temp_update, gud, grid_);

    update_start.ymin++;
  }

  if (update_r.zmin > update_start.zmin)
  {
    temp_update = update_start;
    temp_update.zmax = temp_update.zmin + 1;

    GridUpdateTiling<Update, GridUpdateData, PrivateGridUpdateData>
      ::grid_loop(temp_update, gud, grid_);

    update_start.zmin++;
  }

  if (update_r.xmax < update_start.xmax)
  {
    temp_update = update_start;
    temp_update.xmin = temp_update.xmax - 1;

    GridUpdateTiling<Update, GridUpdateData, PrivateGridUpdateData>
      ::grid_loop(temp_update, gud, grid_);
      
    update_start.xmax--;
  }

  if (update_r.ymax < update_start.ymax)
  {
    temp_update = update_start;
    temp_update.ymin = temp_update.ymax - 1;

    GridUpdateTiling<Update, GridUpdateData, PrivateGridUpdateData>
      ::grid_loop(temp_update, gud, grid_);

    update_start.ymax--;
  }

  if (update_r.zmax < update_start.zmax)
  {
    temp_update = update_start;
    temp_update.zmin = temp_update.zmax - 1;

    GridUpdateTiling<Update, GridUpdateData, PrivateGridUpdateData>
      ::grid_loop(temp_update, gud, grid_);

    update_start.zmax--;
  }

}

void MetaFDTD::update_e(GridUpdateData &gud)
{
  region_t temp_update, update_start;

  switch (mt_)
  {
  case METAFDTD_ONE:
    GridUpdateTiling<ExGridUpdate, GridUpdateData, PrivateGridUpdateData>
      ::grid_loop(grid_->update_ex_r_, gud, *grid_);

    GridUpdateTiling<EyGridUpdate, GridUpdateData, PrivateGridUpdateData>
      ::grid_loop(grid_->update_ey_r_, gud, *grid_);

    GridUpdateTiling<EzGridUpdate, GridUpdateData, PrivateGridUpdateData>
      ::grid_loop(grid_->update_ez_r_, gud, *grid_);

    break;

  case METAFDTD_THREE:
    // Meta programmable field update
    GridUpdateTiling<ElectricGridUpdate, GridUpdateData, PrivateGridUpdateData>
      ::grid_loop(e_update_r_, gud, *grid_);

    component_correct<ExGridUpdate>(e_update_r_, grid_->update_ex_r_, 
                                    gud, *grid_);

    component_correct<EyGridUpdate>(e_update_r_, grid_->update_ey_r_, 
                                    gud, *grid_);

    component_correct<EzGridUpdate>(e_update_r_, grid_->update_ez_r_, 
                                    gud, *grid_);

    break;

  case METAFDTD_SGI_ORIGIN:
    cerr << "SGI Origin specific meta update not implemented!" << endl;
    break;
  }
}

void MetaFDTD::update_h(GridUpdateData &gud)
{
  region_t temp_update;

  switch (mt_)
  {
  case METAFDTD_ONE:
    GridUpdateTiling<HxGridUpdate, GridUpdateData, PrivateGridUpdateData>
      ::grid_loop(grid_->update_hx_r_, gud, *grid_);

    GridUpdateTiling<HyGridUpdate, GridUpdateData, PrivateGridUpdateData>
      ::grid_loop(grid_->update_hy_r_, gud, *grid_);

    GridUpdateTiling<HzGridUpdate, GridUpdateData, PrivateGridUpdateData>
      ::grid_loop(grid_->update_hz_r_, gud, *grid_);

    break;

  case METAFDTD_THREE:
    GridUpdateTiling<MagneticGridUpdate, GridUpdateData, PrivateGridUpdateData>
      ::grid_loop(h_update_r_, gud, *grid_);

    component_correct<HxGridUpdate>(h_update_r_, grid_->update_hx_r_, 
                                    gud, *grid_);

    component_correct<HyGridUpdate>(h_update_r_, grid_->update_hy_r_, 
                                    gud, *grid_);

    component_correct<HzGridUpdate>(h_update_r_, grid_->update_hz_r_, 
                                    gud, *grid_);

    break;

  case METAFDTD_SGI_ORIGIN:
    cerr << "SGI Origin specific meta update not implemented!" << endl;
    break;
  }

}


void MetaFDTD::compute_update_regions()
{
  // The update region that can be visited by all field component
  // update equations.
  e_update_r_ = grid_->update_ex_r_;
  h_update_r_ = grid_->update_hx_r_;

  if (grid_->update_ey_r_.xmin > e_update_r_.xmin)
    e_update_r_.xmin = grid_->update_ey_r_.xmin;

  if (grid_->update_ez_r_.xmin > e_update_r_.xmin)
    e_update_r_.xmin = grid_->update_ez_r_.xmin;

  if (grid_->update_ey_r_.ymin > e_update_r_.ymin)
    e_update_r_.ymin = grid_->update_ey_r_.ymin;

  if (grid_->update_ez_r_.ymin > e_update_r_.ymin)
    e_update_r_.ymin = grid_->update_ez_r_.ymin;

  if (grid_->update_ey_r_.zmin > e_update_r_.zmin)
    e_update_r_.zmin = grid_->update_ey_r_.zmin;

  if (grid_->update_ez_r_.zmin > e_update_r_.zmin)
    e_update_r_.zmin = grid_->update_ez_r_.zmin;

  // E maxs
  if (grid_->update_ey_r_.xmax < e_update_r_.xmax)
    e_update_r_.xmax = grid_->update_ey_r_.xmax;

  if (grid_->update_ez_r_.xmax < e_update_r_.xmax)
    e_update_r_.xmax = grid_->update_ez_r_.xmax;

  if (grid_->update_ey_r_.ymax < e_update_r_.ymax)
    e_update_r_.ymax = grid_->update_ey_r_.ymax;

  if (grid_->update_ez_r_.ymax < e_update_r_.ymax)
    e_update_r_.ymax = grid_->update_ez_r_.ymax;

  if (grid_->update_ey_r_.zmax < e_update_r_.zmax)
    e_update_r_.zmax = grid_->update_ey_r_.zmax;

  if (grid_->update_ez_r_.zmax < e_update_r_.zmax)
    e_update_r_.zmax = grid_->update_ez_r_.zmax;



  // H mins
  if (grid_->update_hy_r_.xmin > h_update_r_.xmin)
    h_update_r_.xmin = grid_->update_hy_r_.xmin;

  if (grid_->update_hz_r_.xmin > h_update_r_.xmin)
    h_update_r_.xmin = grid_->update_hz_r_.xmin;

  if (grid_->update_hy_r_.ymin > h_update_r_.ymin)
    h_update_r_.ymin = grid_->update_hy_r_.ymin;

  if (grid_->update_hz_r_.ymin > h_update_r_.ymin)
    h_update_r_.ymin = grid_->update_hz_r_.ymin;

  if (grid_->update_hy_r_.zmin > h_update_r_.zmin)
    h_update_r_.zmin = grid_->update_hy_r_.zmin;

  if (grid_->update_hz_r_.zmin > h_update_r_.zmin)
    h_update_r_.zmin = grid_->update_hz_r_.zmin;

  // H maxs
  if (grid_->update_hy_r_.xmax < h_update_r_.xmax)
    h_update_r_.xmax = grid_->update_hy_r_.xmax;

  if (grid_->update_hz_r_.xmax < h_update_r_.xmax)
    h_update_r_.xmax = grid_->update_hz_r_.xmax;

  if (grid_->update_hy_r_.ymax < h_update_r_.ymax)
    h_update_r_.ymax = grid_->update_hy_r_.ymax;

  if (grid_->update_hz_r_.ymax < h_update_r_.ymax)
    h_update_r_.ymax = grid_->update_hz_r_.ymax;

  if (grid_->update_hy_r_.zmax < h_update_r_.zmax)
    h_update_r_.zmax = grid_->update_hy_r_.zmax;

  if (grid_->update_hz_r_.zmax < h_update_r_.zmax)
    h_update_r_.zmax = grid_->update_hz_r_.zmax;


#ifdef DEBUG
  cout << "E 3 update region: " << e_update_r_;
  cout << "H 3 update region: " << h_update_r_;

  // This is handled by component correct. 
//   region_t temp_update;
//   temp_update = grid_->update_ex_r_;
//   temp_update.xmax = temp_update.xmin + 1;
//   cout << "Ex 3 correction: " << temp_update;
  
//   temp_update = grid_->update_ey_r_;
//   temp_update.ymax = temp_update.ymin + 1;
//   cout << "Ey 3 correction: " << temp_update;
  
//   temp_update = grid_->update_ez_r_;
//   temp_update.zmax = temp_update.zmin + 1;
//   cout << "Ez 3 correction: " << temp_update;

//   temp_update = grid_->update_hx_r_;
//   temp_update.xmin = temp_update.xmax - 1;
//   cout << "Hx 3 correction: " << temp_update;

//   temp_update = grid_->update_hy_r_;
//   temp_update.ymin = temp_update.ymax - 1;
//   cout << "Hy 3 correction: " << temp_update;

//   temp_update = grid_->update_hz_r_;
//   temp_update.zmin = temp_update.zmax - 1;
//   cout << "Hz 3 correction: " << temp_update;
#endif

}
