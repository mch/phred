/** \class Material 
 * \brief This is a sort of place holder object for material
 * information.
 *
 * It's the go-between from input to grid coefficients. 
 */

#ifndef MATERIAL
#define MATERIAL

#include <string>

using namespace std;

#include "Types.hh"

class Material 
{
 private:
  MaterialType type_;

  /**
   * Electric permittivity 
   */
  mat_prop_t epsilon_;

  /**
   * Conductivity
   */
  mat_prop_t sigma_;
  
  mat_prop_t mu_; /**< Magentic permeability */
  mat_prop_t sigma_star_; /**< The other thing... */

  string name_; /**< Human readable material name */

 protected:

 public:
  
  Material();
  ~Material();
  
  /**
   * Returns the material type
   *
   * @return material type
   */
  MaterialType get_type();

  /**
   * Returns the permittivity
   *
   * @return epsilon
   */
  mat_prop_t get_epsilon();

  /**
   * Returns the conductivity
   *
   * @return sigma
   */
  mat_prop_t get_sigma();

  /**
   * Returns the permeability
   *
   * @return mu
   */
  mat_prop_t get_mu();

  /**
   * Returns the other thing
   *
   * @return sigma star
   */
  mat_prop_t get_sigma_star();

  /**
   * Returns the human readable material name
   *
   * @return name
   */
  string get_name();

  /**
   * Set the material type
   * @param type MaterialType to set
   */
  void set_type(MaterialType type);

  /**
   * Set the permitivity
   * @param eps permitivity
   */
  void set_epsilon(mat_prop_t eps);

  /**
   * Set the conductivity
   * @param sigma conductivity
   */
  void set_sigma(mat_prop_t sigma);

  /**
   * Set the permeability
   * @param mu permeability
   */
  void set_mu(mat_prop_t mu);

  /** 
   * Set the other thing
   * @param sigma_star other thing
   */
  void set_sigma_star(mat_prop_t sigstar);

  /**
   * Set the human readable material name
   * @param name
   */
  void set_name(string name);
};

#endif // MATERIAL
