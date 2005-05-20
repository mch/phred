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

/** \class Material 
 * \brief This is a sort of place holder object for material
 * information.
 *
 * It's the go-between from input to grid coefficients. 
 *
 * The following properties are used to define materials for these
 * dispersions:
 *
 * Debye: 
 *   debye_tau - relaxation time (sec)
 *   debye_eps_inf - permittivity as frequency goes to infinity
 *   debye_eps_s - Permittivity as frequency goes to zero
 *
 * Drude: 
 *   drude_epsilon_inf - Reletive permittivity at infinite frequency
 *   drude_plasma_freq - Plasma frequency (rad/sec)
 *   drude_vc - Collision frequency (Hz)
 *
 * Lorentz: NONE YET
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
  friend class MaterialLib;
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
  bool pec_;

  mat_prop_t mu_; /**< Magentic permeability */
  mat_prop_t sigma_star_; /**< Magnetic loss */

  string name_; /**< Human readable material name */

  map<string, mat_prop_t> properties_; /**< Additional material
                                          properties that may be used
                                          by the grid. */ 

 protected:

 public:
  
  Material();
  ~Material();

  /**
   * Returns a reference to a map of material material properties. 
   */ 
  const map<string, mat_prop_t> &get_properties() const;

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
  void set_property(const char *name, mat_prop_t prop);

  /**
   * Returns a named material property. THROWS AN EXCEPTION if the
   * requested property is not found. 
   */ 
  mat_prop_t get_property(const char *name) const;

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
   * Get the collision frequency
   */
  mat_prop_t get_collision_freq() const;

  /**
   * Get the plasma frequency
   */
  mat_prop_t get_plasma_freq() const;

  /**
   * Determine if this material is PEC
   */ 
  inline bool is_pec() const
  { return pec_; }

  /**
   * Returns the type of material this is by looking at the material
   * properties for clues. For course, PEC and lossy dielectric are
   * easy.
   */ 
  MaterialType type() const;

  /**
   * A conveniece function to help with setting up Drude materials. 
   *
   * @param eps_inf Reletive permittivity at infinite frequency
   * @param plasma_freq Plasma frequency (rad/sec)
   * @param collision_freq Collision frequency (Hz)
   */
  void setup_drude(mat_prop_t eps_inf, mat_prop_t plasma_freq, 
                   mat_prop_t collision_freq);
};

#endif // MATERIAL
