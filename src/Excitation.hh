/** \class Excitation
 * \brief An abstract base class for excitations. 
 *
 * Subclasses contain the data and an implementation of excite(Grid
 * &grid) which knows how to apply the excitation to the grid.
 */

#ifndef EXCITE_H
#define EXCITE_H

#include "Types.hh"
#include "Grid.hh"

class Excitation 
{
protected:
  // Region to apply the source to (in global coordinates)
  unsigned int x_start_;
  unsigned int y_start_;
  unsigned int z_start_;

  unsigned int x_end_;
  unsigned int y_end_;
  unsigned int z_end_;

  /**
   * Field polarization to excite, x, y, z
   */
  field_t polarization_[3];

  /**
   * Field to excite, E or H
   */
  FieldType type_;

  /**
   * True if this is a soft source (defaults to false). Source
   * sources are added to the field rather than the field being set
   * to them.
   */
  bool soft_;

public:
  Excitation();
  virtual ~Excitation();

  /**
   * This function is called to apply this excitation to the grid. 
   *
   * @param grid the grid to operate on
   * @param time_step the time step we are on
   * @param type the field type we are allowed to excite. 
   */
  virtual void excite(Grid &grid, unsigned int time_step, 
                      FieldType type) = 0;

  /**
   * Set the (rectangular) region to apply the excitation to. The
   * start and end coordinates can range from 0 to global grid size -
   * 1. 
   *
   * @param x_start The starting x coordinate
   * @param x_end The ending x coordinate
   * @param y_start The starting y coordinate
   * @param y_end The ending y coordinate
   * @param z_start The starting z coordinate
   * @param z_end The ending z coordinate
   */ 
  void set_region(unsigned int x_start, unsigned int x_end, 
                  unsigned int y_start, unsigned int y_end, 
                  unsigned int z_start, unsigned int z_end);

  /**
   * Set the polarization vector
   *
   * @param x the x component
   * @param y the y component
   * @param z the z component
   */
  void set_polarization(field_t x, field_t y, field_t z);

  /**
   * Set the field to excite. Either E or H. 
   *
   * @param t FieldType
   */
  void set_type(FieldType t);

  /**
   * Tells if this source is "soft" or not. Soft sources are
   * additive. 
   */
  inline bool get_soft()
  {
    return soft_;
  }

  /**
   * Set the softness of the source.
   * @param soft true if this source should be soft.
   */
  inline void set_soft(bool soft)
  {
    soft_ = soft;
  }

};

#endif // EXCITE_H

