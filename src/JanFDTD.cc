#include "JanFDTD.hh"

/**
 * A function implemented in the yacc file. 
 */ 
void parse_jan_grammer(const char *filename, JanFDTD *jfdtd);

JanFDTD::JanFDTD()
{}

JanFDTD::~JanFDTD()
{
  map<string, Excitation *>::iterator iter;
  map<string, Excitation *>::iterator iter_e = e_excitations_.end();

  for (iter = e_excitations_.begin(); iter != iter_e; ++iter)
    delete iter->second;

  iter_e = h_excitations_.end();
  for (iter = h_excitations_.begin(); iter != iter_e; ++iter)
    delete iter->second;

  map<string, Result *>::iterator riter;
  map<string, Result *>::iterator riter_e = results_.end();

  for(riter = results_.begin(); riter != riter_e; ++riter)
    delete riter->second;

  map<string, DataWriter *>::iterator diter;
  map<string, DataWriter *>::iterator diter_e = datawriters_.end();

  for(diter = datawriters_.begin(); diter != diter_e; ++diter)
    delete diter->second;

  vector<Geometry *>::iterator giter;
  vector<Geometry *>::iterator giter_e = geometry_.end();

  for(giter = geometry_.begin(); giter != giter_e; ++giter)
    delete *giter;

  if (mlib_)
    delete mlib_;
}

void JanFDTD::parse_file(string filename)
{
  parse_jan_grammer(filename.c_str(), this);
}

void JanFDTD::set_program_mode(const char *m)
{

}

void JanFDTD::set_structure_mode(const char *m)
{

}

void JanFDTD::set_timestep_mode(const char *m)
{

}

void JanFDTD::set_runtime(unsigned int rt)
{

}

void JanFDTD::set_time_modulo(unsigned int rt)
{

}

void JanFDTD::set_dimx(unsigned int dx)
{

}

void JanFDTD::set_dimy(unsigned int dy)
{

}

void JanFDTD::set_dimz(unsigned int dz)
{

}


void JanFDTD::set_deltax(double dx)
{

}

void JanFDTD::set_deltay(double dy)
{

}

void JanFDTD::set_deltaz(double dz)
{

}

void JanFDTD::set_deltat(double dt)
{

}


void JanFDTD::set_num_materials(unsigned int nm)
{

}

void JanFDTD::add_material(unsigned int id, float eps, float sigma, 
                  float mu, float sigma_star)
{

}


void JanFDTD::set_pml_boundary(Face face, unsigned int thickness, 
                      char p, unsigned int nrml_refl)
{

}


void JanFDTD::set_boundary(Face face, BoundaryCondition bc)
{

}

