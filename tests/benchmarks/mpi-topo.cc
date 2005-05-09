// Benchmarks to explore the difference that using MPI topologies to
// map processes to the hardware makes. One would expect that mapping
// processes to the hardware such that adjacent processes are
// physically close to each other will reduce communication time.

#include <mpi.h>

#include <iostream>
#include <fstream>

using namespace std;

static const int LOOPS = 1000;
static int MIN_CUBE_SZ = 20;
static int MAX_CUBE_SZ = 300;
static const int INCR = 1;
static bool use_topo_g = false;
static char *hostname_g;

static int
decode_switches (int argc, char **argv)
{
  opterr = 0;
  int c;
  char *arg = 0;

  while ((c = getopt (argc, argv, "a:b:xyz")) != -1)
    switch (c)
    {
    case 'a':
      arg = const_cast<char *>(optarg);

      MIN_CUBE_SZ = atoi(arg);
      break;

    case 'b':
      arg = const_cast<char *>(optarg);
      
      MAX_CUBE_SZ = atoi(arg);
      break;

    case 'x':
      run_yz_g = true;
      break;

    case 'y':
      run_zx_g = true;
      break;

    case 'z':
      run_xy_g = true;
      break;
    }

  return 0;
}

static void
usage()
{
  cout << "mpi-topo: MPI topology benchmarks for FDTD style communication\n\
\n\
Usage: mpi-topo [-a MIN_CUBE_SIZE] [-b MAX_CUBE_SIZE] [-t]\n\
\n\
 -a Set the minimum cube size to use. Defaults to 20.\n\
 -b Set the maximum cube size to use. Defaults to 200.\n\
 -t Use MPI tolopology to map processes to hardware (defaults to random).\n\
\n\n";

}

typedef struct {
  double bandwidth;
  double time;
  double latency;
  double bytes;
  double message_size; 
} result_t;

void set_comm()
{

}

int main(int argc, char **argv)
{


  return 0;
}
