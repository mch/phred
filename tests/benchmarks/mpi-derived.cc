/**
 * A simple bandwidth (bytes per second) and latency (sec) measurment
 * tool for non-contiguous data of various sizes transmitted using MPI
 * derived datatypes. 
 *
 * The data is a 3d cubic space, and the data transmitted is XY, YZ,
 * ZX planes of data within that cube. This mimics the transactions
 * that are used in Phred.
 */ 

#include <mpi.h>
#include <iostream>
#include <fstream>

#include <unistd.h>
#include <stdlib.h>

using namespace std;

static const int LOOPS = 1000;
static int MIN_CUBE_SZ = 20;
static int MAX_CUBE_SZ = 300;
static const int INCR = 1;
static bool run_yz_g = false;
static bool run_zx_g = false;
static bool run_xy_g = false;
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
  cout << "mpi-derived: Simple MPI benchmarks using derived datatypes\n\
\n\
Usage: mpi-derived [-a MIN_CUBE_SIZE] [-b MAX_CUBE_SIZE] [-x] [-y] [-z]\n\
\n\
 -a Set the minimum cube size to use. Defaults to 20.\n\
 -b Set the maximum cube size to use. Defaults to 200.\n\
 -x Test MPI transmission in the X=Constant (yz) plane.\n\
 -y Test MPI transmission in the Y=Constant (zx) plane.\n\
 -z Test MPI transmission in the Z=Constant (xy) plane.\n\
\n\
At least one of -x, -y, or -z must be selected. \n\n";

}

typedef struct {
  double bandwidth;
  double time;
  double latency;
  double bytes;
  double message_size; 
} result_t;

result_t bandwidth(double *buf, MPI::Datatype &dt)
{
  double start, stop;
  result_t r;
  MPI::Status status;

  if (MPI::COMM_WORLD.Get_rank() == 0)
  {
    start = MPI::Wtime();
    for (int i = 0; i < LOOPS; i++)
    {
      MPI::COMM_WORLD.Send(buf, 1, dt, 1, 1);
      MPI::COMM_WORLD.Recv(buf, 1, dt, 1, 1, status);
    }
    stop = MPI::Wtime();
  } 
  else 
  {
    start = MPI::Wtime();
    for (int i = 0; i < LOOPS; i++)
    {
      MPI::COMM_WORLD.Recv(buf, 1, dt, 0, 1, status);
      MPI::COMM_WORLD.Send(buf, 1, dt, 0, 1);
    }
    stop = MPI::Wtime();
  }

  r.time = stop - start;
  r.bytes = 2 * dt.Get_size() * LOOPS;
  r.bandwidth = r.bytes / r.time;
  r.latency = r.time / LOOPS;
  r.message_size = dt.Get_size();

  return r;
}

void loop_xy(double *buf)
{
  ofstream bwfile;
  MPI::Datatype xy_plane;
  MPI::Datatype y_vector;
  result_t r;
  string fn(hostname_g);
  fn += "_xy_bw.txt";

  if (MPI::COMM_WORLD.Get_rank() == 0)
    bwfile.open(fn.c_str(), ofstream::out);

  for (int sz = MIN_CUBE_SZ; sz <= MAX_CUBE_SZ; sz += INCR)
  {
    y_vector = MPI::DOUBLE.Create_vector(sz, 1, sz);

    xy_plane = y_vector.Create_hvector(sz, 1, sizeof(double) * sz * sz);

    xy_plane.Commit();

    r = bandwidth(buf, xy_plane);

    if (MPI::COMM_WORLD.Get_rank() == 0)
      bwfile << sz << " " << r.time << " " << r.bandwidth 
             << " " << r.latency << " " << r.bytes 
             << " " << r.message_size << endl;

    xy_plane.Free();
  }

  if (MPI::COMM_WORLD.Get_rank() == 0)
    bwfile.close();

}

void loop_yz(double *buf)
{
  ofstream bwfile;
  MPI::Datatype yz_plane;
  result_t r;
  string fn(hostname_g);
  fn += "_yz_bw.txt";

  if (MPI::COMM_WORLD.Get_rank() == 0)
    bwfile.open(fn.c_str(), ofstream::out);

  for (int sz = MIN_CUBE_SZ; sz <= MAX_CUBE_SZ; sz += INCR)
  {
    yz_plane = MPI::DOUBLE.Create_contiguous(sz * sz);

    yz_plane.Commit();

    r = bandwidth(buf, yz_plane);

    if (MPI::COMM_WORLD.Get_rank() == 0)
      bwfile << sz << " " << r.time << " " << r.bandwidth 
             << " " << r.latency << " " << r.bytes 
             << " " << r.message_size << endl;

    yz_plane.Free();
  }

  if (MPI::COMM_WORLD.Get_rank() == 0)
    bwfile.close();

}

void loop_zx(double *buf)
{
  ofstream bwfile;
  MPI::Datatype zx_plane;
  result_t r;
  string fn(hostname_g);
  fn += "_zx_bw.txt";
  
  if (MPI::COMM_WORLD.Get_rank() == 0)
    bwfile.open(fn.c_str(), ofstream::out);

  for (int sz = MIN_CUBE_SZ; sz <= MAX_CUBE_SZ; sz += INCR)
  {
    zx_plane = MPI::DOUBLE.Create_vector(sz, sz, sz * sz);

    zx_plane.Commit();

    r = bandwidth(buf, zx_plane);

    if (MPI::COMM_WORLD.Get_rank() == 0)
      bwfile << sz << " " << r.time << " " << r.bandwidth 
             << " " << r.latency << " " << r.bytes 
             << " " << r.message_size << endl;

    zx_plane.Free();
  }

  if (MPI::COMM_WORLD.Get_rank() == 0)
    bwfile.close();

}

int main(int argc, char **argv)
{
  double *buf;

  MPI::Datatype yz_plane;
  MPI::Datatype zx_plane;

  MPI::Init(argc, argv);

  decode_switches(argc, argv);

  hostname_g = new char[255];
  gethostname(hostname_g, 255);

  if (run_yz_g || run_zx_g || run_xy_g)
  {

    if (MPI::COMM_WORLD.Get_size() == 2)
    {
      buf = new double[MAX_CUBE_SZ * MAX_CUBE_SZ * MAX_CUBE_SZ];
    
      if (run_yz_g)
        loop_yz(buf);
      
      if (run_xy_g)
        loop_xy(buf);

      if (run_zx_g)
        loop_zx(buf);

      delete[] buf;
    } else {
      cout << "This benchmark only works with two mpi processes. " << endl;
      usage();
    }
  } else {
    usage();
  }

  delete[] hostname_g;

  MPI::Finalize();
  
  return 0;
}
