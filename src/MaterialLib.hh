/** \class MaterialLib
 * \brief A material library.
 * This class is a material library for keeping track of materials. 
 * It contains a vector of material objects, and is able to generate
 * arrays of material coefficients for use in the FDTD algorithm. 
 */

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

  /**
   * Add a material to the material library. A copy of the material
   * object is stored, not a pointer or a reference.
   *
   * @param mat the material to store
   */
  void add_material(Material& mat);

  /**
   * Returns the number of materials in the library
   * @return an int
   */
  int num_materials();

  /**
   * Returns an interator that points to the start of the vector
   * containing the materials. TEMPORARY until a proper iterator for
   * the materials is set up.
   *
   * @return a vector<Material>::iterator
   */
  vector<Material>::iterator get_material_iter_begin();

  /**
   * Returns an interator that points to the end of the vector
   * containing the materials. TEMPORARY until a proper iterator for
   * the materials is set up.
   *
   * @return a vector<Material>::iterator
   */
  vector<Material>::iterator get_material_iter_end();
};

#endif // MATERIAL_LIB
