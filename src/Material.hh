// This is a sort of place holder object for material
// information. It's the go-between from input to grid coefficients. 

#ifndef MATERIAL
#define MATERIAL

#include <string>

using namespace std;

#include "Types.hh"

class Material 
{
 private:
  MaterialType type_;

  // Electric permittivity and conductivtiy
  mat_prop_t epsilon_;
  mat_prop_t sigma_;
  
  // Magentic permeability and whatever the other magnetic thing is. 
  mat_prop_t mu_;
  mat_prop_t sigma_star_;

  // Human readable material name
  string name_;

 protected:

 public:
  
  Material();
  ~Material();
  
  MaterialType get_type();
  mat_prop_t get_epsilon();
  mat_prop_t get_sigma();
  mat_prop_t get_mu();
  mat_prop_t get_sigma_star();
  string get_name();

  void set_type(MaterialType type);
  void set_epsilon(mat_prop_t eps);
  void set_sigma(mat_prop_t sigma);
  void set_mu(mat_prop_t mu);
  void set_sigma_star(mat_prop_t sigstar);
  void set_name(string name);
};

#endif // MATERIAL
