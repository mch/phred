
// This class is a material library for keeping track of materials. 
// It contains a vector of material objects, and is able to generate
// arrays of material coefficients for use in the FDTD algorithm. 

#ifndef MATERIAL_LIB
#define MATERIAL_LIB

#include <vector>
#include <iostream>
#include <string>

#include "Types.hh"
#include "Material.hh"

using namespace std; 

class MaterialLib {
 private:
  vector <Material> materials_;

 protected:

 public:
  // An input iterator class for traversing the materials list
  //class iterator {
  //private:
  //  vector<Material>::iterator iter_;
  //  
  //public:
  //  
  //};

  // Functions
  MaterialLib();
  ~MaterialLib();

  void add_material(Material& mat);

  int num_materials();

  // Temp, until I finish the iterator above
  vector<Material>::iterator get_material_iter_begin();
  vector<Material>::iterator get_material_iter_end();
};

#endif // MATERIAL_LIB
