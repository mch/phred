/**
 * A simple bandwidth (bytes per second) and latency (sec) measurment
 * tool for contiguous data of various sizes.
 */ 

#include <mpi.h>
#include <iostream>
#include <fstream>

#include <unistd.h>
#include <stdlib.h>

using namespace std;

static const int LOOPS = 1000;
static int MIN_BYTES = 1000;
static int MAX_BYTES = 200000;
static const int INCR = 1000;
static char *hostname_g;
static bool simple_g = false;
static bool pingpong_g = false;

static int
decode_switches (int argc, char **argv)
{
  opterr = 0;
  int c;
  char *arg = 0;

  while ((c = getopt (argc, argv, "a:b:sp")) != -1)
    switch (c)
    {
    case 'a':
      arg = const_cast<char *>(optarg);

      MIN_BYTES = atoi(arg);
      break;

    case 'b':
      arg = const_cast<char *>(optarg);
      
      MAX_BYTES = atoi(arg);
      break;

    case 's':
      simple_g = true;
      break;

    case 'p':
      pingpong_g = true;
      break;
    }

  return 0;
}

static void
usage()
{
  cout << "mpi-simple: Simple MPI benchmarks using contiguous data\n\
\n\
Usage: mpi-simple [-a MIN_BYTES] [-b MAX_BYTES] [-s] [-p]\n\
\n\
 -a Set the minimum message size (bytes). Defaults to 1000.\n\
 -b Set the maximum message size (bytes). Defaults to 200000.\n\
 -s Test simple one-way MPI transmission.\n\
 -p Test ping-pong MPI transmission.\n\
\n\
At least one of -s, or -p must be selected. \n\n";

}

typedef struct {
  double bandwidth;
  double time;
  double latency;
} result_t;

result_t simple(char *buf, unsigned int message_size)
{
  double start, stop;
  result_t r;

  if (MPI::COMM_WORLD.Get_rank() == 0)
  {
    start = MPI::Wtime();
    for (int i = 0; i < LOOPS; i++)
    {
      MPI::COMM_WORLD.Send(buf, message_size, MPI::CHAR, 1, 1);
    }
    stop = MPI::Wtime();
  } 
  else 
  {
    MPI::Status status;
    start = MPI::Wtime();
    for (int i = 0; i < LOOPS; i++)
    {
      MPI::COMM_WORLD.Recv(buf, message_size, MPI::CHAR, 0, 1, status);
    }
    stop = MPI::Wtime();
  }

  r.time = stop - start;
  r.bandwidth = (message_size * LOOPS) / r.time;
  r.latency = r.time / LOOPS;

  return r;
}

result_t pingpong(char *buf, unsigned int message_size)
{
  double start, stop;
  result_t r;
  MPI::Status status;

  if (MPI::COMM_WORLD.Get_rank() == 0)
  {
    start = MPI::Wtime();
    for (int i = 0; i < LOOPS; i++)
    {
      MPI::COMM_WORLD.Send(buf, message_size, MPI::CHAR, 1, 1);
      MPI::COMM_WORLD.Recv(buf, message_size, MPI::CHAR, 1, 1, status);
    }
    stop = MPI::Wtime();
  } 
  else 
  {
    start = MPI::Wtime();
    for (int i = 0; i < LOOPS; i++)
    {
      MPI::COMM_WORLD.Recv(buf, message_size, MPI::CHAR, 0, 1, status);
      MPI::COMM_WORLD.Send(buf, message_size, MPI::CHAR, 0, 1);
    }
    stop = MPI::Wtime();
  }

  r.time = stop - start;
  r.bandwidth = (2 * message_size * LOOPS) / r.time;
  r.latency = r.time / LOOPS;

  return r;
}

int main(int argc, char **argv)
{
  char *buf;
  
  MPI::Init(argc, argv);

  hostname_g = new char[255];
  gethostname(hostname_g, 255);

  decode_switches(argc, argv);

  if (MPI::COMM_WORLD.Get_size() == 2)
  {
    buf = new char[MAX_BYTES];

    if (simple_g)
    {
      ofstream bwfile;
      string fn(hostname_g);
      fn += "_simple.txt";
      
      if (MPI::COMM_WORLD.Get_rank() == 0)
        bwfile.open(fn.c_str(), ofstream::out);
      
      result_t r;
      for (int sz = MIN_BYTES; sz <= MAX_BYTES; sz += INCR)
      {
        r = simple(buf, sz);
        
        if (MPI::COMM_WORLD.Get_rank() == 0)
          bwfile << sz << " " << r.time << " " << r.bandwidth 
                 << " " << r.latency << endl;
      }
      
      if (MPI::COMM_WORLD.Get_rank() == 0)
        bwfile.close();
    }

    // ugg, cut and paste
    if (pingpong_g)
    {
      ofstream bwfile;
      string fn(hostname_g);
      fn += "_pingpong.txt";
      
      if (MPI::COMM_WORLD.Get_rank() == 0)
        bwfile.open(fn.c_str(), ofstream::out);
      
      result_t r;
      for (int sz = MIN_BYTES; sz <= MAX_BYTES; sz += INCR)
      {
        r = simple(buf, sz);
        
        if (MPI::COMM_WORLD.Get_rank() == 0)
          bwfile << sz << " " << r.time << " " << r.bandwidth 
                 << " " << r.latency << endl;
      }
      
      if (MPI::COMM_WORLD.Get_rank() == 0)
        bwfile.close();
    }

    delete[] buf;

    if (!simple_g && !pingpong_g)
    {
      cout << "At least one of -s or -p must be selected. " << endl;
      usage();
    }

  } else {
    cout << "This benchmark only works with two mpi processes. " << endl;
    usage();
  }

  delete[] hostname_g;

  MPI::Finalize();
  
  return 0;
}
