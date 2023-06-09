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

#include <boost/python.hpp>

#include "../Grid.hh"
#include "../Boundaries/Ewall.hh"
#include "../Boundaries/Hwall.hh"
#include "../Boundaries/UPml.hh"
#include "../Boundaries/Pml.hh"
#include "../Boundaries/Periodic.hh"
#include "../Types.hh"

using namespace boost::python;

void call_apply(BoundaryCond &bc, Face face, Grid &grid, FieldType type)
{
  return bc.apply(face, grid, type);
}

BoundaryCondition call_get_type(BoundaryCond &bc)
{
  return bc.get_type();
}

void call_init(BoundaryCond &bc, const Grid &grid, Face face)
{
  return bc.init(grid, face);
}

void call_deinit(BoundaryCond &bc, const Grid &grid, Face face)
{
  return bc.deinit(grid, face);
}

void call_add_sd_bcs(BoundaryCond &bc, SubdomainBc *sd, 
                     Face bcface, Face sdface)
{
  return bc.add_sd_bcs(sd, bcface, sdface);
}

/**
 * Boundary condition wrapper so that boundary conditions can be
 * defined in Python.
 */ 
class BoundaryCondWrap : public BoundaryCond
{
private:
  PyObject *self_;

public:
  BoundaryCondWrap(PyObject *self)
    : self_(self)
  {}
  
  void apply(Face face, Grid &grid, FieldType type)
  {
    return call_method<void>(self_, "apply");
  }

  BoundaryCondition get_type() const
  {
    return call_method<BoundaryCondition>(self_, "get_type");
  }

  void init(const Grid &grid, Face face)
  {
    return call_method<void>(self_, "init");
  }

  void deinit(const Grid &grid, Face face)
  {
    return call_method<void>(self_, "deinit");
  }

  void add_sd_bcs(SubdomainBc *sd, Face bcface, Face sdface)
  {
    return call_method<void>(self_, "add_sd_bcs");
  }

};

void export_boundaries()
{
  class_<BoundaryCond, BoundaryCondWrap, boost::noncopyable>("BoundaryCond", "Boundary Condition base class")
    .add_property("thickness", &BoundaryCond::get_thickness, 
                  &BoundaryCond::set_thickness 
                  /*, "Number of cells thick the boundary is." */)
    .def("apply", &BoundaryCondWrap::apply, 
         "Apply the boundary condition. Should probably never be called.")
    .def("get_type", &BoundaryCondWrap::get_type, 
         "Returns an integer representing the type of boundary this is.")
    .def("init", &BoundaryCondWrap::init, 
         "Initialize the boundary. Should probably never be called.")
    .def("deinit", &BoundaryCondWrap::deinit, 
         "De-initialize the boundary. Should probably never be called.")
    ;

  class_<Ewall, bases<BoundaryCond> >("Ewall", "Electric wall")
    .def("apply", &Ewall::apply)
    .def("get_type", &Ewall::get_type)
    ;

  class_<Hwall, bases<BoundaryCond> >("Hwall", "Magnetic wall")
    .def("apply", &Hwall::apply)
    .def("get_type", &Hwall::get_type)
    ;

  class_<Pml, bases<BoundaryCond> >("Pml", "Berenger's split field PML")
    .def("apply", &Pml::apply)
    .def("set_variation", &Pml::set_variation)
    .def("set_g_param", &Pml::set_g_param)
    .def("set_nrml_relf", &Pml::set_nrml_refl)
    .def("get_type", &Pml::get_type)
    .def("set_thickness", &UPml::set_thickness,
         "Number of cells thick the PML is.")
    ;

  class_<UPml, bases<BoundaryCond>, boost::noncopyable >
    ("UPml", "Gedney's uniaxial PML")
    .def("apply", &UPml::apply)
    .def("get_type", &UPml::get_type)
    .def("set_thickness", &UPml::set_thickness,
         "The number of cells occupied by the PML.")
    .add_property("sigma_max", &UPml::get_sigma_max, 
                  &UPml::set_sigma_max //, 
                  //"Maximum conductivity. Calculated by default."
                  )
    .add_property("poly_order", &UPml::get_poly_order, 
                  &UPml::set_poly_order //, 
                  //"Order of the polynomial grading. Defaults to 4."
                  )
    .add_property("eps_opt", &UPml::get_eps_opt,
                  &UPml::set_eps_opt //, 
                  //"Relative permittivity for the purpose of calculating "
                  //"the maximum conductivity of the PML. If there are big "
                  //"permittivity continuities withing the PML, the best "
                  //"performace may be to average them. Defaults to 1.0." 
                  )
    .add_property("sigma_ratio", &UPml::get_sigma_ratio, 
                  &UPml::set_sigma_ratio //,
                  //"The ratio between sigma_opt, calculated from the "
                  //"polynomial order and material properties, to "
                  //"sigma_max, the maximum condictivity of the UPML. "
                  //"Defaults to 1.0."
                  )
    .add_property("k_max", &UPml::get_k_max,
                  &UPml::set_k_max //,
                  //"This parameter can be set to a value > 1 to "
                  //"terminate evanescent waves. It is graded into the "
                  //"PML like sigma max. Setting this to 30 had been "
                  //"found useful."
                  )
    ;

  class_<Periodic, bases<BoundaryCond>, boost::noncopyable>
    ("Periodic", "Periodic boundary condition",
     init<shared_ptr<PeriodicExcitation> >())
    .def("apply", &Periodic::apply)
    ;
                                     
}
