// This file defines the types used for various things. This
// may be processed by autoconf in the future so these things 
// can be configure switches, or dependent on machine type. 

#ifndef TYPES_H
#define TYPES_H

/**
 * Boundary conditions. These can be appled to faces of the grid or to
 * the surfaces of materials. But normally just to the faces of the
 * grid.
 */
enum BoundaryCondition {
  UNKNOWN, /**< The boundary condition is unknown; non-sane. */
  SUBDOMAIN, /**< Means that the face is shared with another node */
  EWALL, /**< Electric wall */
  HWALL, /**< Magnetic wall */
  PML /**< Perfectly matched layers (Berenger's absorbing boundary) */
  //MUR /**< 2nd order absorbing boundary condition */
  //IMPEDANCE, /**< Imedance boundary */
};

/** 
 * Material types. Not currently used, but may be used to indicate the
 * type of grid that should be used.
 */
enum MaterialType {
  PERF_COND, 
  NON_PERMEABLE, 
  CONDUCTIVE, 
  DEBYE, 
  LORENTZ, 
  DRUDE
};

/**
 * Faces on the grid cube.
 */
enum Face {
  FRONT = 0, /**< x = dimx, YZ plane */
  BACK = 1, /**< x = 0, YZ plane */
  LEFT = 2, /**< y = 0, XZ plane */
  RIGHT = 3, /**< y = dimy, XZ plane */
  BOTTOM = 4, /**< z = 0, XY plane */
  TOP = 5 /**< z = dimz, XY plane */
};

/**
 * Field components. Electric and magnetic
 */
enum FieldComponent {
  EX,
  EY,
  EZ,
  HX, 
  HY, 
  HZ
};

/**
 * Field types, E and H
 */
enum FieldType {
  E,
  H
};

/**
 * The type used for material coefficients
 */
typedef float mat_coef_t; 

/**
 * The type used for material properties
 */
typedef double mat_prop_t;

/**
 * The type used for grid field components, THE MPI DATATYPE MUST MATCH!
 */
typedef double field_t;
#define GRID_MPI_TYPE MPI_DOUBLE

/**
 * The type used for grid spacings
 */
typedef double delta_t;

/** \struct region
 * \brief min and max coordinates that can be used to define a
 * block within the grid. Or a face. 
 */
typedef struct region_t {
  unsigned int xmin;
  unsigned int xmax;
  unsigned int ymin;
  unsigned int ymax;
  unsigned int zmin;
  unsigned int zmax;
} region;

#endif // TYPES_H
