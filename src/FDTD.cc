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


#include "FDTD.hh"

// Globals from phred.cc
bool mnps, estimate_memory;

FDTD::FDTD()
  : grid_(0), mlib_(0)
{}

FDTD::~FDTD()
{
  if (grid_)
    delete grid_;
}

Geometry &FDTD::find_geometry(unsigned int x,
                              unsigned int y,
                              unsigned int z)
{
  vector<Geometry *>::reverse_iterator iter;
  vector<Geometry *>::reverse_iterator iter_e = geometry_.rend();

  for (iter = geometry_.rbegin(); iter != iter_e; ++iter)
  {
    if ((*iter)->local_point_inside(x, y, z))
      return **iter;
  }
}

void FDTD::set_time_steps(unsigned int t)
{
  time_steps_ = t;
}

void FDTD::set_grid_size(unsigned int x, unsigned int y,
                         unsigned int z)
{
  global_ginfo_.global_dimx_ = global_ginfo_.dimx_ = x;
  global_ginfo_.global_dimy_ = global_ginfo_.dimy_ = y;
  global_ginfo_.global_dimz_ = global_ginfo_.dimz_ = z;
}

void FDTD::set_grid_deltas(field_t dx, field_t dy, field_t dz)
{
  global_ginfo_.deltax_ = dx;
  global_ginfo_.deltay_ = dy;
  global_ginfo_.deltaz_ = dz;
  
  global_ginfo_.deltat_ = 0.9 / 
    ( C * sqrt( 1/(pow(dx, static_cast<float>(2.0))) + 
                1/(pow(dy, static_cast<float>(2.0))) + 
                1/(pow(dz, static_cast<float>(2.0)))));
}

void FDTD::set_time_delta(field_t dt)
{
  global_ginfo_.deltat_ = dt;
}

void FDTD::set_boundary(Face face, BoundaryCond *bc)
{
  global_ginfo_.set_boundary(face, bc);
}

void FDTD::load_materials(MaterialLib &matlib)
{
  mlib_ = &matlib;
}

void FDTD::add_excitation(const char *name, Excitation *ex)
{
  FieldType t = ex->get_type();
  if (t == E || t == BOTH)
    e_excitations_[string(name)] = ex;

  if (t == H || t == BOTH)
    h_excitations_[string(name)] = ex;    
}

void FDTD::add_result(const char *name, Result *r)
{
  results_[string(name)] = r;
  //if (r->get_name().length() == 0)
  r->set_name(name);
}

void FDTD::add_datawriter(const char *name, DataWriter *dw)
{
  datawriters_[string(name)] = dw;
}

void FDTD::add_geometry(Geometry *g)
{
  geometry_.push_back(g);
}

void FDTD::map_result_to_datawriter(const char *result, const char *dw)
{
  map<string, Result *>::iterator riter = results_.find(result);
  map<string, DataWriter *>::iterator dwiter = datawriters_.find(dw);

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
    map<string, Result *>::iterator riter 
      = results_.find((*iter).first);

    map<string, DataWriter *>::iterator dwiter
      = datawriters_.find((*iter).second);

    if (riter != results_.end() && dwiter != datawriters_.end())
    {
      (*dwiter).second->add_variable(*(*riter).second);
    }
    
    ++iter;
  }
}

// CHOP THIS UP; make helper functions
void FDTD::run(int rank, int size)
{
  // Grid setup
  SimpleSDAlg dd;

  local_ginfo_ = dd.decompose_domain(rank, size, global_ginfo_);

  // Decide what grid to used from materials
  bool freqgrid = false;
  vector<Material>::const_iterator miter = mlib_->get_material_iter_begin();
  vector<Material>::const_iterator miter_e = mlib_->get_material_iter_end();
  
  while(miter != miter_e)
  {
    if ((*miter).get_collision_freq() != 0.0 ||
        (*miter).get_plasma_freq() != 0.0)
    {
      freqgrid = true; 
      break;
    }
    miter++;
  }

  if (grid_)
    delete grid_;

  if (freqgrid)
    {
      grid_ = new FreqGrid();
      cout << "Using freq grid. " << endl;
    }
  else
    {
      cout << "Using simple grid. " << endl;
      grid_ = new Grid();
    }

  grid_->setup_grid(local_ginfo_);

  grid_->load_materials(*mlib_);
  grid_->load_geometries(geometry_);

  grid_->set_define_mode(false);

  // Life cycle init
  map<string, Excitation *>::iterator e_eiter_b = e_excitations_.begin();
  map<string, Excitation *>::iterator e_eiter = e_eiter_b;
  map<string, Excitation *>::iterator e_eiter_e = e_excitations_.end();
  
  while (e_eiter != e_eiter_e) 
  {
    (*e_eiter).second->init(*grid_);
    ++e_eiter;
  }

  map<string, Excitation *>::iterator h_eiter_b = h_excitations_.begin();
  map<string, Excitation *>::iterator h_eiter = h_eiter_b;
  map<string, Excitation *>::iterator h_eiter_e = h_excitations_.end();
  
  while (h_eiter != h_eiter_e) 
  {
    (*h_eiter).second->init(*grid_);
    ++h_eiter;
  }

  map<string, Result *>::iterator riter_b = results_.begin();
  map<string, Result *>::iterator riter = riter_b;
  map<string, Result *>::iterator riter_e = results_.end();
  
  while (riter != riter_e) 
  {
    (*riter).second->init(*grid_);
    ++riter;
  }

  map<string, DataWriter *>::iterator dwiter_b = datawriters_.begin();
  map<string, DataWriter *>::iterator dwiter = dwiter_b;
  map<string, DataWriter *>::iterator dwiter_e = datawriters_.end();
  
  while (dwiter != dwiter_e) 
  {
    (*dwiter).second->init(*grid_);
    ++dwiter;
  }

  vector<Geometry *>::iterator giter = geometry_.begin();
  vector<Geometry *>::iterator giter_e = geometry_.end();
  
  while(giter != giter_e)
  {
    (*giter)->init(*grid_);
    (*giter)->set_material(*grid_);
    ++giter;
  }
  
  setup_datawriters();

  // Run
  unsigned int ts = 0;

  // For optionally tracking millions of nodes per second. 
  time_t start, now, time_total = 0;
  clock_t start_cpu, now_cpu, time_total_cpu = 0;

  for (ts = 1; ts <= time_steps_; ts++) {
    cout << "phred time step " << ts << endl;
    
    // Fields update
    if (mnps) 
    {
      start=time(NULL);
      start_cpu = clock();
    }

    grid_->update_h_field();

    if (mnps) 
    {
      now = time(NULL);
      now_cpu = clock();
      time_total += now - start;
      time_total_cpu += now_cpu - start_cpu;
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

    // Fields update
    if (mnps) 
    {
      start=time(NULL);
      start_cpu = clock();
    }

    grid_->update_e_field();
    
    if (mnps) 
    {
      now = time(NULL);
      now_cpu = clock();
      time_total += now - start;
      time_total_cpu += now_cpu - start_cpu;
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
    vector< pair<string, string> >::iterator iter = r_dw_map_.begin();
    vector< pair<string, string> >::iterator iter_e = r_dw_map_.end();

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

  if (mnps)
  {
    time_t avg_time = time_total / time_steps_;
    clock_t avg_cpu_time = time_total_cpu / time_steps_;
    unsigned int num_mnodes = grid_->get_num_updated_nodes() / 1e6;

    cout << "Average millions of nodes per second, w.r.t. wall clock time: " 
         << num_mnodes / avg_time << endl;
    cout << "Average millions of nodes per second, w.r.t. CPU time: " 
         << num_mnodes / (avg_cpu_time / static_cast<double>(CLOCKS_PER_SEC)) 
         << endl;
  }

  // life cycle de init
  e_eiter = e_eiter_b;
  while (e_eiter != e_eiter_e) 
  {
    (*e_eiter).second->deinit(*grid_);
    ++e_eiter;
  }

  h_eiter = h_eiter_b;
  while (h_eiter != h_eiter_e) 
  {
    (*h_eiter).second->deinit(*grid_);
    ++h_eiter;
  }

  riter = riter_b;
  while (riter != riter_e) 
  {
    (*riter).second->deinit(*grid_);
    ++riter;
  }

  dwiter = dwiter_b;
  while (dwiter != dwiter_e) 
  {
    (*dwiter).second->deinit(*grid_);
    ++dwiter;
  }

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

