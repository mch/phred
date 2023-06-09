/**
 * A simple bandwidth (bytes per second) and latency (sec) measurment
 * tool for non-contiguous data of various sizes transmitted using
 * manual buffer packing and unpacking.
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
  cout << "mpi-manual: Simple MPI benchmarks using manual packing\n\
\n\
Usage: mpi-manual [-a MIN_CUBE_SIZE] [-b MAX_CUBE_SIZE] [-x] [-y] [-z]\n\
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

void loop_xy(double *buf)
{
  ofstream bwfile;
  result_t r;
  string fn(hostname_g);
  fn += "_xy_bw.txt";

  if (MPI::COMM_WORLD.Get_rank() == 0)
    bwfile.open(fn.c_str(), ofstream::out);

  for (int sz = MIN_CUBE_SZ; sz <= MAX_CUBE_SZ; sz += INCR)
  {
    int msgsz = sz * sz; 
    int stride = sz;
    double *temp = new double[sz * sz];
    int buf_offset = 0;
    int temp_offset = 0;

    double start, stop;
    result_t r;
    MPI::Status status;

    if (MPI::COMM_WORLD.Get_rank() == 0)
    {
      start = MPI::Wtime();
      for (int i = 0; i < LOOPS; i++)
      {
        // Pack 
        for (int xi = 0; xi < msgsz; xi++)
        {
          temp[temp_offset] = buf[buf_offset];
          temp_offset++;
          buf_offset += stride;
        }

        // send/recv
        MPI::COMM_WORLD.Send(temp, msgsz, MPI::DOUBLE, 1, 1);
        MPI::COMM_WORLD.Recv(temp, msgsz, MPI::DOUBLE, 1, 1, status);

        // Unpack
        buf_offset = 0;
        temp_offset = 0;
        for (int xi = 0; xi < msgsz; xi++)
        {
          buf[buf_offset] = temp[temp_offset];
          temp_offset++;
          buf_offset += stride;
        }
      }
      stop = MPI::Wtime();
    } 
    else 
    {
      start = MPI::Wtime();
      for (int i = 0; i < LOOPS; i++)
      {
        MPI::COMM_WORLD.Recv(temp, msgsz, MPI::DOUBLE, 0, 1, status);

        // Unpack
        buf_offset = 0;
        temp_offset = 0;
        for (int xi = 0; xi < msgsz; xi++)
        {
          buf[buf_offset] = temp[temp_offset];
          temp_offset++;
          buf_offset += stride;
        }

        // Pack 
        buf_offset = 0;
        temp_offset = 0;
        for (int xi = 0; xi < msgsz; xi++)
        {
          temp[temp_offset] = buf[buf_offset];
          temp_offset++;
          buf_offset += stride;
        }

        MPI::COMM_WORLD.Send(temp, msgsz, MPI::DOUBLE, 0, 1);
      }
      stop = MPI::Wtime();
    }

    r.time = stop - start;
    r.bytes = 2 * msgsz * LOOPS;
    r.bandwidth = r.bytes / r.time;
    r.latency = r.time / LOOPS;
    r.message_size = msgsz;
    
    if (MPI::COMM_WORLD.Get_rank() == 0)
      bwfile << sz << " " << r.time << " " << r.bandwidth 
             << " " << r.latency << " " << r.bytes 
             << " " << r.message_size << endl;

    delete[] temp;
  }

  if (MPI::COMM_WORLD.Get_rank() == 0)
    bwfile.close();

}

// No packing needed
void loop_yz(double *buf)
{
  ofstream bwfile;
  result_t r;
  string fn(hostname_g);
  fn += "_yz_bw.txt";

  if (MPI::COMM_WORLD.Get_rank() == 0)
    bwfile.open(fn.c_str(), ofstream::out);

  for (int sz = MIN_CUBE_SZ; sz <= MAX_CUBE_SZ; sz += INCR)
  {
    int msgsz = sz * sz; 
    double start, stop;
    result_t r;
    MPI::Status status;

    if (MPI::COMM_WORLD.Get_rank() == 0)
    {
      start = MPI::Wtime();
      for (int i = 0; i < LOOPS; i++)
      {
        MPI::COMM_WORLD.Send(buf, msgsz, MPI::DOUBLE, 1, 1);
        MPI::COMM_WORLD.Recv(buf, msgsz, MPI::DOUBLE, 1, 1, status);
      }
      stop = MPI::Wtime();
    } 
    else 
    {
      start = MPI::Wtime();
      for (int i = 0; i < LOOPS; i++)
      {
        MPI::COMM_WORLD.Recv(buf, msgsz, MPI::DOUBLE, 0, 1, status);
        MPI::COMM_WORLD.Send(buf, msgsz, MPI::DOUBLE, 0, 1);
      }
      stop = MPI::Wtime();
    }

    r.time = stop - start;
    r.bytes = 2 * msgsz * LOOPS;
    r.bandwidth = r.bytes / r.time;
    r.latency = r.time / LOOPS;
    r.message_size = msgsz;
    
    if (MPI::COMM_WORLD.Get_rank() == 0)
      bwfile << sz << " " << r.time << " " << r.bandwidth 
             << " " << r.latency << " " << r.bytes 
             << " " << r.message_size << endl;
  }

  if (MPI::COMM_WORLD.Get_rank() == 0)
    bwfile.close();

}

void loop_zx(double *buf)
{
  ofstream bwfile;
  result_t r;
  string fn(hostname_g);
  fn += "_zx_bw.txt";
  
  if (MPI::COMM_WORLD.Get_rank() == 0)
    bwfile.open(fn.c_str(), ofstream::out);

  for (int sz = MIN_CUBE_SZ; sz <= MAX_CUBE_SZ; sz += INCR)
  {
    int msgsz = sz * sz; 
    int stride = sz * sz;
    double *temp = new double[sz * sz];

    double start, stop;
    result_t r;
    MPI::Status status;

    if (MPI::COMM_WORLD.Get_rank() == 0)
    {
      start = MPI::Wtime();
      for (int i = 0; i < LOOPS; i++)
      {
        // Pack

        // Send / Recv
        MPI::COMM_WORLD.Send(temp, msgsz, MPI::DOUBLE, 1, 1);
        MPI::COMM_WORLD.Recv(temp, msgsz, MPI::DOUBLE, 1, 1, status);

        // Unpack

        
      }
      stop = MPI::Wtime();
    }
    else 
    {
      start = MPI::Wtime();
      for (int i = 0; i < LOOPS; i++)
      {
        MPI::COMM_WORLD.Recv(temp, msgsz, MPI::DOUBLE, 0, 1, status);
        MPI::COMM_WORLD.Send(temp, msgsz, MPI::DOUBLE, 0, 1);
      }
      stop = MPI::Wtime();
    }

    r.time = stop - start;
    r.bytes = 2 * msgsz * LOOPS;
    r.bandwidth = r.bytes / r.time;
    r.latency = r.time / LOOPS;
    r.message_size = msgsz;
    
    if (MPI::COMM_WORLD.Get_rank() == 0)
      bwfile << sz << " " << r.time << " " << r.bandwidth 
             << " " << r.latency << " " << r.bytes 
             << " " << r.message_size << endl;

    delete[] temp;
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
