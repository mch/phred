// This file defines the types used for various things. This
// may be processed by autoconf in the future so these things 
// can be configure switches, or dependent on machine type. 

#ifndef TYPES_H
#define TYPES_H

// Boundary conditions
enum BoundaryCondition {
  SUBDOMAIN, // Share me
  EWALL, 
  HWALL, 
  ESYM,
  HSYM,
  IMPEDANCE,
  PML
};

// Material types
enum MaterialType {
  PERF_COND, 
  NON_PERMEABLE, 
  CONDUCTIVE, 
  DEBYE, 
  LORENTZ, 
  DRUDE
};

// Faces
enum Face {
  FRONT = 0, 
  BACK = 1, 
  LEFT = 2, 
  RIGHT = 3,
  BOTTOM = 4, 
  TOP = 5
};

// Field components
enum FieldComponent {
  EX,
  EY,
  EZ,
  HX, 
  HY, 
  HZ
};

// Material coefficients
typedef float mat_coef_t; 

// Material properties
typedef double mat_prop_t;

// Grid field components, THE MPI DATATYPE MUST MATCH!
typedef double field_t;
#define GRID_MPI_TYPE MPI_DOUBLE

// Grid spacings
typedef double delta_t;


#endif // TYPES_H
