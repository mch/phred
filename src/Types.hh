// This file defines the types used for various things. This
// may be processed by autoconf in the future so these things 
// can be configure switches, or dependent on machine type. 

// Material types
enum MaterialType {PERF_COND}

// Material coefficients
typedef float mat_coef_t; 

// Material properties
typedef double mat_prop_t;

// Grid field components, THE MPI DATATYPE MUST MATCH!
typedef double field_t;
#define GRID_MPI_TYPE MPI_DOUBLE
