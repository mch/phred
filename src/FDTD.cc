
#include "FDTD.hh"

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
  
  global_ginfo_.deltat_ = 1 / 
    ( C * sqrt( 1/(pow(dx, static_cast<float>(2.0))) + 
                1/(pow(dy, static_cast<float>(2.0))) + 
                1/(pow(dz, static_cast<float>(2.0)))));
}

void FDTD::set_boundary(Face face, BoundaryCond *bc)
{
  global_ginfo_.set_boundary(face, bc);
}

void FDTD::load_materials(MaterialLib &matlib)
{
  grid_.load_materials(matlib);
}

void FDTD::add_excitation(const char *name, Excitation *ex)
{
  excitations_[string(name)] = ex;
}

void FDTD::add_result(const char *name, Result *r)
{
  results_[string(name)] = r;
}

void FDTD::add_datawriter(const char *name, DataWriter *dw)
{
  datawriters_[string(name)] = dw;
}

void FDTD::map_result_to_datawriter(const char *result, const char *dw)
{
  map<string, Result *>::iterator riter = results_.find(result);
  map<string, DataWriter *>::iterator dwiter = datawriters_.find(dw);

  if (riter != results_.end() && dwiter != datawriters_.end())
    r_dw_map_.push_back(pair<string, string>(string(result), string(dw)));
  else
    throw FDTDException("Result and data writer must be present before mapping them together");
}

void FDTD::run(unsigned int steps)
{
  // Grid setup

  // Life cycle init
  map<string, Excitation *>::iterator eiter_b = excitations_.begin();
  map<string, Excitation *>::iterator eiter = eiter_b;
  map<string, Excitation *>::iterator eiter_e = excitations_.end();
  
  while (eiter != eiter_e) 
  {
    (*eiter).second->init(grid_);
    ++eiter;
  }

  map<string, Result *>::iterator riter_b = results_.begin();
  map<string, Result *>::iterator riter = riter_b;
  map<string, Result *>::iterator riter_e = results_.end();
  
  while (riter != riter_e) 
  {
    (*riter).second->init(grid_);
    ++riter;
  }

  map<string, DataWriter *>::iterator dwiter_b = datawriters_.begin();
  map<string, DataWriter *>::iterator dwiter = dwiter_b;
  map<string, DataWriter *>::iterator dwiter_e = datawriters_.end();
  
  while (dwiter != dwiter_e) 
  {
    (*dwiter).second->init(grid_);
    ++dwiter;
  }

  // Run
  unsigned int ts = 0;

  for (ts = 1; ts < steps; ts++) {
    cout << "phred time step " << ts << endl;

    // Fields update
    grid_.update_h_field();

    // Boundary condition application
    grid_.apply_boundaries(H);

    // Excitations
    eiter = eiter_b;
    while (eiter != eiter_e)
    {
      (*eiter).second->excite(grid_, ts, H);
      ++eiter;
    }

    // Fields update
    grid_.update_e_field();
    
    // Boundary condition application
    grid_.apply_boundaries(E);

    // Excitations
    eiter = eiter_b;
    while (eiter != eiter_e)
    {
      (*eiter).second->excite(grid_, ts, E);
      ++eiter;
    }
    
    // Results
    vector< pair<string, string> >::iterator iter = r_dw_map_.begin();
    vector< pair<string, string> >::iterator iter_e = r_dw_map_.end();

    while (iter != iter_e)
    {
      riter = results_.find((*iter).first);
      dwiter = datawriters_.find((*iter).second);      
      
      if (riter != riter_e && dwiter != dwiter_e)
        (*dwiter).second->handle_data(ts, 
                                      (*riter).second->get_result(grid_, ts));

      ++iter;
    }
  }

  // life cycle de init
  eiter = eiter_b;
  while (eiter != eiter_e) 
  {
    (*eiter).second->deinit(grid_);
    ++eiter;
  }

  riter = riter_b;
  while (riter != riter_e) 
  {
    (*riter).second->deinit(grid_);
    ++riter;
  }

  dwiter = dwiter_b;
  while (dwiter != dwiter_e) 
  {
    (*dwiter).second->deinit(grid_);
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

