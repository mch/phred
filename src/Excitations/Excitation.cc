/* 
   Phred - Phred is a parallel finite difference time domain
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

#include "Excitation.hh"

Excitation::Excitation(shared_ptr<Signal> sf)
  : type_(E), soft_(false),
    sf_(sf), time_start_(0.0), time_stop_(0.0), time_offset_(0.0)
{
  polarization_[0] = 1;
  polarization_[1] = 0;
  polarization_[2] = 0;
}

Excitation::~Excitation()
{}

void Excitation::set_region(shared_ptr<CSGBox> box)
{
  box_ = box;
}

void Excitation::set_polarization(field_t x, field_t y, field_t z)
{
  polarization_[0] = x;
  polarization_[1] = y;
  polarization_[2] = z;
}

void Excitation::set_type(FieldType t)
{
  type_ = t;
}

void Excitation::init(const Grid &grid)
{
  if (box_.get())
  {
    region_ = grid.get_local_region(*(box_.get()));

#ifdef DEBUG
    cout << "Excitation local region: (" << (*region_).xmin() << ", "
         << (*region_).ymin() << ", " << (*region_).zmin() << ") -> ("
         << (*region_).xmax() << ", " << (*region_).ymax() << ", " 
         << (*region_).zmax() << ")." << endl;
#endif
  } else {
    throw RegionException("Excitations must be assigned a region "
                          "over which they are to be applied!");
  }
}

void Excitation::excite(Grid &grid, unsigned int time_step, 
                        FieldType type)
{
  // Find out where we fit in this grid (convert to local coordinates)
  if (type != BOTH && type != type_)
    return;

  field_t e_time = grid.get_deltat() * time_step;
  field_t h_time = grid.get_deltat() * (time_step - 0.5);
  field_t e_sf = sf_->signal_function(e_time - time_offset_);
  field_t h_sf = sf_->signal_function(h_time - time_offset_);

  if (h_time < time_start_ || (e_time > time_stop_ 
                               && time_stop_ > grid.get_deltat()))
    return;

  field_t e_fld[3];
  field_t h_fld[3];

  e_fld[0] = e_sf * polarization_[0];
  e_fld[1] = e_sf * polarization_[1];
  e_fld[2] = e_sf * polarization_[2];

  h_fld[0] = h_sf * polarization_[0];
  h_fld[1] = h_sf * polarization_[1];
  h_fld[2] = h_sf * polarization_[2];

  if (!soft_) 
  {
    for(unsigned int i = (*region_).xmin(); i < (*region_).xmax(); i++)
    {
      for (unsigned int j = (*region_).ymin(); j < (*region_).ymax(); j++)
      {
        for (unsigned int k = (*region_).zmin(); k < (*region_).zmax(); k++)
        {
          switch (type_) 
          {
          case E:
            grid.set_ex(i,j,k, e_fld[0]);
            grid.set_ey(i,j,k, e_fld[1]);
            grid.set_ez(i,j,k, e_fld[2]);
            break;

          case H:
            grid.set_hx(i,j,k, h_fld[0]);
            grid.set_hy(i,j,k, h_fld[1]);
            grid.set_hz(i,j,k, h_fld[2]);
            break;

          case BOTH: // Isn't meant for Excitations.
            throw std::exception();
            break; 
          }
        }
      }
    }
  } else {
    for(unsigned int i = (*region_).xmin(); i < (*region_).xmax(); i++)
    {
      for (unsigned int j = (*region_).ymin(); j < (*region_).ymax(); j++)
      {
        for (unsigned int k = (*region_).zmin(); k < (*region_).zmax(); k++)
        {
          switch (type_) 
          {
          case E:
            if (polarization_[0] != 0.0) grid.set_ex(i,j,k, 
                                                     e_fld[0] + grid.get_ex(i,j,k));
            if (polarization_[1] != 0.0) grid.set_ey(i,j,k, 
                                                     e_fld[1] + grid.get_ey(i,j,k));
            if (polarization_[2] != 0.0) grid.set_ez(i,j,k, 
                                                     e_fld[2] + grid.get_ez(i,j,k));
            break;

          case H:
            if (polarization_[0] != 0.0) grid.set_hx(i,j,k, 
                                                     h_fld[0] + grid.get_hx(i,j,k));
            if (polarization_[1] != 0.0) grid.set_hy(i,j,k, 
                                                     h_fld[1] + grid.get_hy(i,j,k));
            if (polarization_[2] != 0.0) grid.set_hz(i,j,k, 
                                                     h_fld[2] + grid.get_hz(i,j,k));
            break;

          case BOTH: // Isn't meant for Excitations.
            throw std::exception();
            break; 
          }
        }
      }
    }
  }
}

ostream& Excitation::to_string(ostream &os) const
{
  os << "In the local grid, this excitation is being applied to "
     << (*region_) << ". The polarization vector is ("
     << polarization_[0] << ", " << polarization_[1]
     << ", " << polarization_[2] << "). The excitation is ";

  if (soft_)
    os << "soft";
  else
    os << "hard";

  os << " and is being applied to the ";

  if (type_ == E)
    os << "electric field. ";
  else if (type_ == H)
    os << "magnetic field. ";
  else 
    os << "electric and magnetic fields.";

  os << "The signal function being applied is " << (*sf_);

  return os;
}
