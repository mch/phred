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

/** \class Material 
 * \brief This is a sort of place holder object for material
 * information.
 *
 * It's the go-between from input to grid coefficients. 
 */

#ifndef MATERIAL
#define MATERIAL

#include <string>
#include <map>

using namespace std;

#include "Types.hh"
#include "Exceptions.hh"

class Material 
{
 private:
  MaterialType type_;
  unsigned int material_id_; /**< id used to reference constants in
                                the grid. */

  /**
   * Electric permittivity 
   */
  mat_prop_t epsilon_;

  /**
   * Conductivity
   */
  mat_prop_t sigma_;
  
  mat_prop_t mu_; /**< Magentic permeability */
  mat_prop_t sigma_star_; /**< Magnetic loss */

  string name_; /**< Human readable material name */

  // These are used for plasma implementation
  mat_prop_t vc_; /**< Collision frequency */
  mat_prop_t fp_; /**< Plasma frequency */

  map<string, mat_prop_t> properties_; /**< Additional material
                                          properties that may be used
                                          by the grid. */ 

 protected:

 public:
  
  Material();
  ~Material();
  
  /**
   * Returns the material ID used to index into arrays of
   * constants. Grid use only. 
   */ 
  inline unsigned int get_id() const
  {
    return material_id_;
  }

  /**
   * Used by the Grid to set the id of this material. 
   */ 
  inline void set_id(unsigned int id)
  {
    material_id_ = id;
  }

  /**
   * Set a named material property. Any existing value is
   * overwritten. Consult the documentation for the Grid classes for
   * a list of material properties required. 
   */ 
  inline void set_property(const char *name, mat_prop_t prop)
  {
    properties_[name] = prop;
  }

  /**
   * Returns a named material property. THROWS AN EXCEPTION if the
   * requested property is not found. 
   */ 
  inline mat_prop_t get_property(const char *name) const
  {
    map<string, mat_prop_t>::const_iterator iter = properties_.find(name);
    
    if (iter == properties_.end())
      throw MaterialPropertyException("Property not found.");

    return (*iter).second;
  }

  /**
   * Returns the material type
   *
   * @return material type
   */
  MaterialType get_type() const;

  /**
   * Returns the permittivity
   *
   * @return epsilon
   */
  mat_prop_t get_epsilon() const;

  /**
   * Returns the conductivity
   *
   * @return sigma
   */
  mat_prop_t get_sigma() const;

  /**
   * Returns the permeability
   *
   * @return mu
   */
  mat_prop_t get_mu() const;

  /**
   * Returns the other thing
   *
   * @return sigma star
   */
  mat_prop_t get_sigma_star() const;

  /**
   * Returns the human readable material name
   *
   * @return name
   */
  string get_name() const;

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
  void set_name(const char *name);

  /**
   * Set the collision frequency for plasma materials
   */
  void set_collision_freq(mat_prop_t vc);

  /**
   * Set the plasma frequency, in radians per second!
   */
  void set_plasma_freq(mat_prop_t fp);
  
  /**
   * Get the collision frequency
   */
  mat_prop_t get_collision_freq() const;

  /**
   * Get the plasma frequency
   */
  mat_prop_t get_plasma_freq() const;

};

#endif // MATERIAL
