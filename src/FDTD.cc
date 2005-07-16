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


#include "FDTD.hh"
#include "Globals.hh"

FDTD::FDTD()
  : time_steps_(0), dt_scale_(0.9)
{
  mlib_ = shared_ptr<MaterialLib>(new MaterialLib()); // Empty default. 
  dd_alg_ = shared_ptr<SubdomainAlg>(new MPISubdomainAlg());
}

FDTD::~FDTD()
{}

void FDTD::set_grid_material(const char *material)
{
  geometry_.set_grid_material(material);
}

unsigned int FDTD::get_num_x_cells()
{
  return global_ginfo_.global_dimx_;
}

unsigned int FDTD::get_num_y_cells()
{
  return global_ginfo_.global_dimy_;
}

unsigned int FDTD::get_num_z_cells()
{
  return global_ginfo_.global_dimz_;
}

void FDTD::set_time_steps(unsigned int t)
{
  time_steps_ = t;
}

void FDTD::set_time(float t)
{
  time_steps_ = static_cast<unsigned int>(t / global_ginfo_.deltat_);
}

void FDTD::set_grid_size(float x, float y, float z)
{
  geometry_.set_grid_size(x, y, z);

  if (global_ginfo_.deltax_ > 0)
    global_ginfo_.global_dimx_ = global_ginfo_.dimx_ = 
      static_cast<unsigned int>(floor(x / global_ginfo_.deltax_));

  if (global_ginfo_.deltay_ > 0)
    global_ginfo_.global_dimy_ = global_ginfo_.dimy_ = 
      static_cast<unsigned int>(floor(y / global_ginfo_.deltay_));

  if (global_ginfo_.deltaz_ > 0)
    global_ginfo_.global_dimz_ = global_ginfo_.dimz_ = 
      static_cast<unsigned int>(floor(z / global_ginfo_.deltaz_));
}

void FDTD::set_grid_centre(float x, float y, float z)
{
  geometry_.set_grid_centre(x, y, z);
}

void FDTD::set_grid_deltas(field_t dx, field_t dy, field_t dz)
{
  global_ginfo_.deltax_ = dx;
  global_ginfo_.deltay_ = dy;
  global_ginfo_.deltaz_ = dz;
  
  global_ginfo_.deltat_ = dt_scale_ / 
    ( C * sqrt( 1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz)));

  point gsize = geometry_.get_grid_size();

  if (global_ginfo_.deltax_ > 0)
    global_ginfo_.global_dimx_ = global_ginfo_.dimx_ = 
      static_cast<unsigned int>(ceil(gsize.x / global_ginfo_.deltax_));

  if (global_ginfo_.deltay_ > 0)
    global_ginfo_.global_dimy_ = global_ginfo_.dimy_ = 
      static_cast<unsigned int>(ceil(gsize.y / global_ginfo_.deltay_));

  if (global_ginfo_.deltaz_ > 0)
    global_ginfo_.global_dimz_ = global_ginfo_.dimz_ = 
      static_cast<unsigned int>(ceil(gsize.z / global_ginfo_.deltaz_));
}

field_t FDTD::get_time_delta()
{
  return global_ginfo_.deltat_;
}

void FDTD::set_time_delta(field_t dt)
{
  global_ginfo_.deltat_ = dt;
}

void FDTD::set_boundary(Face face, shared_ptr<BoundaryCond> bc)
{
  global_ginfo_.set_boundary(face, bc);
}

void FDTD::load_materials(shared_ptr<MaterialLib> matlib)
{
  mlib_ = matlib;
}

void FDTD::add_excitation(const char *name, shared_ptr<Excitation> ex)
{
  FieldType t = ex->get_type();
  if (t == E || t == BOTH)
    e_excitations_[string(name)] = ex;

  if (t == H || t == BOTH)
    h_excitations_[string(name)] = ex;    
}

void FDTD::add_result(const char *name, shared_ptr<Result> r)
{
  results_[string(name)] = r;
  //if (r->get_name().length() == 0)
  r->set_name(name);
}

void FDTD::add_result(const char *name, shared_ptr<Result> r,
                      const char *dw)
{
  add_result(name, r);
  map_result_to_datawriter(name, dw);
}

void FDTD::add_datawriter(const char *name, shared_ptr<DataWriter> dw)
{
  datawriters_[string(name)] = dw;
}

void FDTD::add_object(string material, shared_ptr<CSGObject> obj)
{
  geometry_.add_object(material, obj);
}

void FDTD::map_result_to_datawriter(const char *result, const char *dw)
{
  map<string, shared_ptr<Result> >::iterator 
    riter = results_.find(result);
  map<string, shared_ptr<DataWriter> >::iterator 
    dwiter = datawriters_.find(dw);

  if (riter != results_.end() && dwiter != datawriters_.end())
  {
    r_dw_map_.push_back(pair<string, string>(string(result), string(dw)));
  }
  else
    throw FDTDException("Result and data writer must be present before mapping them together");
}

void FDTD::setup_datawriters()
{
  vector< pair<string, string> >::iterator iter = r_dw_map_.begin();  
  vector< pair<string, string> >::iterator iter_e = r_dw_map_.end();  

  while (iter != iter_e)
  {
    map<string, shared_ptr<Result> >::iterator riter 
      = results_.find((*iter).first);

    map<string, shared_ptr<DataWriter> >::iterator dwiter
      = datawriters_.find((*iter).second);

    if (riter != results_.end() && dwiter != datawriters_.end())
    {
      (*dwiter).second->add_variable(*(riter->second));
    }
    
    ++iter;
  }
}

// CHOP THIS UP; make helper functions
void FDTD::run()
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

  if (dd_alg_.get())
    local_ginfo_ = dd_alg_->decompose_domain(global_ginfo_);
  else {
    cout << "No domain decomposition algorithm available! Aborting.\n";
    return;
  }

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
  {
    // It's nice to know what's going on
    grid_->get_grid_info().to_string(cout);

    cout << "Initializing grid..." << endl;
  }

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
    
      if (MPI_RANK == 0 && (ts - 10) % 100 == 0)
      {
        rt_now = time(NULL);
        int secs = static_cast<int>((static_cast<double>(rt_now - rt_start) 
                                     / rt_steps) * (time_steps_ - ts));

        cout << "Estimated time remaining at time step " << ts << ": \n\t";
        print_elapsed_time(secs);

        rt_steps = 100;
        rt_start = time(NULL);
      }

      grid_->update_h_field();

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

      grid_->update_e_field();

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

void FDTD::print_elapsed_time(int secs)
{
  int mins = secs / 60;
  secs = secs % 60;
  int hours = mins / 60;
  mins = mins % 60;
  int days = hours / 24;
  hours = hours % 24;
  
  if (days == 1)
    cout << days << " day, "; 
  else if (days > 1)
    cout << days << " days, "; 
  
  if (hours == 1)
    cout << hours << " hour, ";
  else if (hours > 1)
    cout << hours << " hours, ";
  
  if (mins == 1)
    cout << mins << " minute, ";
  else if (mins > 1)
    cout << mins << " minutes, ";
  
  cout << secs << " seconds. " << endl;
}


// template<class T, class A>
// void FDTD::init_objs()
// {
//   map<string, T>::iterator iter = A.begin();
//   map<string, T>::iterator iter_e = A.end();
  
//   while (iter != iter_e) 
//   {
//     (*iter).second->init(grid_);
//     ++iter;
//   }
// }

// template<class T, class A>
// void FDTD::deinit_objs()
// {
//   map<string, T>::iterator iter = A.begin();
//   map<string, T>::iterator iter_e = A.end();
  
//   while (iter != iter_e) 
//   {
//     (*iter).second->deinit(grid_);
//     ++iter;
//   }

// }

