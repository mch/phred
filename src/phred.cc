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


/** \mainpage Phred Documentation
 *
 * \section Introduction
 *
 * Phred is a parallel finite difference time domain electromagnetics
 * simulator. It is written in C and C++, and uses MPI (Message
 * passing interface) for parallelism. It runs best on a cluster of
 * homogenous nodes, where the number of nodes is an even number. 
 *
 * \section addtocode Adding to the code
 * 
 * Phred is designed to be an extensible and flexible research
 * code. As such, it is hopefully easy for other programmers to add
 * boundary conditions, excitations, output formats, and other cool
 * things. This section summarizes how to add such items to Phred in
 * C++. Phred embeds the Python language, and such items may also be
 * prototyped in Python by creating derived classes, although for
 * speed the implementation should be done in C/C++. 
 *
 * \subsection addbc Adding a boundary condition
 *
 * Boundary conditions can be added by subclassing the
 * BoundaryCondition object. For boundary conditions which only need
 * to consider the tangential and normal field components at a face,
 * adapter classes can be used to simplify the expression of the
 * boundary condition. See the electric (Ewall) and magnetic wall
 * boundary conditions (Hwall) for details. 
 *
 * \subsection adde Adding an excitation
 *
 * Phred has built in support for excitations which vary only with
 * time and for excitations which vary with time and
 * space. Excitations, by default, are defined over some rectangular
 * region of space. It's possible to change this by overriding the
 * Excite::excite() method of the Excitation class. 
 *
 * The classes TimeExcitation and SpaceExcitation take care of looping
 * over the default block of space by implementing the excite()
 * method. Sources which need to be applied to some other shape of
 * space can override Excitation::excite() to do so. 
 *
 * Any new excitation which is content to be applied to a block of the
 * grid needs to only override the source_function() method of
 * TimeExcitation or SpaceExcitation. 
 *
 * The source_function() method is public, so that other things, such
 * as the total/scattered thing can use the same excitation objects. 
 *
 *
 * \subsection addoutput Adding a output method
 *
 * Phred supports MATLAB friendly ASCII and binary, netCDF, and HDF
 * output formats. To add support for you favorite file format, make a
 * class which inherits from ... and implement the methods of the
 * interface. 
 */

/* MPI (rocks your socks right off) */
#include <mpi.h>

/* Let's use C++ for things that aren't speed critical, because life
   is just so much easier that way. And safer. Practice safe hex. */
#include <iostream>
#include <vector>
#include <string>

#ifdef USE_RUSAGE
/* rusage() */ 
#include <sys/ctime>
#include <sys/resource.h>
#endif

using namespace std; // Too lazy to type namespaces all the time. 

#include "config.h"

#ifdef USE_PY_BINDINGS
#include <Python.h>
#include "PythonBindings/PyInterpreter.hh"
#endif

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <boost/lexical_cast.hpp>
using boost::lexical_cast;
using boost::bad_lexical_cast;

/* For handling low memory conditions, which happens quite often */
#include <new>

void no_memory()
{
  //throw MemoryException();
  cerr << "Out of memory. Aborting program.";
  MPI_Abort(MPI_COMM_WORLD, 1);
}

// TESTING ONLY! REMOVE AT SOME POINT
#include "Tests.hh"

static void usage (int status);

#ifdef HAVE_LIBPOPT

/* popt plays way nicer with MPI than getopt. Trust me. */
#include <popt.h>

static struct poptOption options[] = 
  {
    {"help", 'h', POPT_ARG_NONE, 0, 'h'},
    {"version", 'V', POPT_ARG_NONE, 0, 'V'},
    {"interactive", 'i', POPT_ARG_NONE, 0, 'i'},
    {"file", 'f', POPT_ARG_STRING, 0, 'f'},
    {"memory", 'm', POPT_ARG_NONE, 0, 'm'},
    {"mnps", 'b', POPT_ARG_NONE, 0, 'b'}, 
    {"quiet", 'q', POPT_ARG_NONE, 0, 'q'}, 
    {"test", 't', POPT_ARG_STRING, 0, 't'},
    {0, 0, 0, 0, 0}
  };
#endif

static int decode_switches (int argc, char **argv);

static string get_extension(string filename);

// Ugly globals
string inputfile, test_case;
const char *program_name;
bool interactive, estimate_memory, mnps, quiet, test_run;
int MPI_RANK, MPI_SIZE, argi_g;

/* Set all the option flags according to the switches specified.
   Return the index of the first non-option argument.  */
static int
decode_switches (int argc, char **argv)
{
#ifdef HAVE_LIBPOPT
  int c;
  char *arg = 0;
  argi_g = 0;

  poptContext ctx = poptGetContext(0, argc, 
                                   const_cast<const char **>(argv), 
                                   options, 0);

  while ((c = poptGetNextOpt (ctx)) > 0 || c == POPT_ERROR_BADOPT)
  {
    argi_g++;

    if (c == POPT_ERROR_BADOPT)
      continue;

    switch (c)
    {
    case 'V':
      cout << "phred " << VERSION << endl;
      exit (0);
      break;

    case 'h':
      usage (0);
      break;

    case 'i':
      interactive = true;
      break;
      
    case 'f':
      arg = const_cast<char *>(poptGetOptArg(ctx));

      if (arg)
        inputfile = arg;
      else 
      {
        cout << "No filename given for the -f switch." << endl;
        usage(0);
      }
      break;

    case 'm':
      estimate_memory = true; 
      break;

    case 'b':
      mnps = true;
      break;

    case 'q':
      quiet = true;
      break;

    case 't':
      arg = const_cast<char *>(poptGetOptArg(ctx));

      if (arg)
        test_case = arg;
      else 
      {
        cout << "No testcase specified." << endl;
        usage(0);
      }

      test_run = true;
      break;

    default:
      cout << "WARNING: got unknown option number: " << c << endl;
    }
  }

  poptFreeContext(ctx);
#else
  mnps = true;
  test_run = true;
#endif

  return 0;
}

static void
usage (int status)
{
  printf ("%s - \
Phred is a parallel finite difference time domain electromagnetics simulator.\n", program_name);

#ifdef HAVE_LIBPOPT
  printf ("Usage: %s [OPTION]... [FILE]...\n", program_name);
  printf ("\
Options:\n\
  -f, --filename             filename to read problem description from\n");
#ifdef USE_PY_BINDINGS
  printf("  -i, --interactive          start an interactive Python interpreter on\n\
                             rank zero if that process is attached to a \n\
                             terminal. Commands will be mirrored to\n\
                             interpreters running on the other ranks.\n");
#endif
  printf("  -h, --help                 display this help and exit\n\
  -q, --quiet                Don't echo configuration information at\n\
                             start up, don't report each time step, etc.\n\
  -m, --memory               Estimate amount of required memory and exit\n\
  -b, --mnps                 Benchmark: estimate the millions of nodes \n\
                             processed per second (this does NOT include\n\
                             time spent in IO or other activities; only\n\
                             node update times are counted.\n\
  -t, --test                 Run a hard coded test problem; select from:\n\
                             H   Single circular hole\n\
                             M   Million node benchmark\n\
                             V   Variable number of nodes benchmark; \n\
                                 follow by the number of cells\n\
                                 in the X, Y, and Z axis.\n\
                             S   Square hole, followed by the size of \n\
                                 the hole in the y dimension as\n\
                                 an integer number of nanometers.\n\
                             through a single hole in a PEC plate.\n\
  -V, --version              output version information and exit\n\
");
#else
  printf ("Usage: %s A [opts]\n", program_name);
  printf("A is a letter identifying the test to run, and [opts] are any\n"
         "options that test takes.\nTests:\n\
  H   Single circular hole\n\
  M   Million node benchmark\n\
  V   Variable number of nodes benchmark; follow by the number of cells\n\
      in the X, Y, and Z axis.\n\
  S   Square hole, followed by the size of the hole in the y dimension as\n\
      an integer number of nanometers.\n\
");
#endif

  exit(status);
}

string get_extension(string filename)
{
  int pos = filename.rfind('.');
  string ext;

  if (pos != string::npos)
    ext = filename.substr(pos + 1, filename.length() - pos);

  return ext;
}

// MAIN!
int main (int argc, char **argv)
{
  time_t start, now;
  double time_total_cpu = 0.0;
  
  // Install a handler for low memory conditions. 
  std::set_new_handler(no_memory);

  start=time(NULL);

  string prog_name;

  interactive = false;
  estimate_memory = false;
  mnps = false;
  quiet = false;

  //std::set_terminate (__gnu_cxx::__verbose_terminate_handler);

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);
  MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);

  //if (MPI_RANK == 0)
  {
    prog_name = argv[0];
    program_name = prog_name.c_str();
    
    // MPI implementations are not required to distribute command line
    // args, although MPICH does.
    decode_switches (argc, argv);
  } 
  // else { // rank 0 passes appropriate args
  
  cout << PACKAGE_NAME << " version " << PACKAGE_VERSION 
       << ", Copyright (C) 2004 Matt Hughes <mch@ieee.org>\n"
       << PACKAGE_NAME << " comes with ABSOLUTELY NO WARRANTY.\n"
       << "This is free software, and you are welcome to redistribute it\n"
       << "under certian conditions. See the COPYING file for details.\n";

  if (interactive)
  {
    cout << "For warranty information type `print warranty'.\n"
         << "For redistribution conditions type `print conditions'.\n";
  }

  cout << "\nMPI information: \nThis process is rank number " << MPI_RANK
       << ".\nThere are a total of " << MPI_SIZE
       << " processes in this group." << endl;

#ifdef USE_OPENMP
  cout << "\nOpenMP information: \nNumber of threads in team: " 
       << omp_get_num_threads()
       << "\nMaximum number of threads in team: " << omp_get_max_threads()
       << "\nNumber of processors: " << omp_get_num_procs()
       << "\nCurrent thread number: " << omp_get_thread_num()
       << "\nDynamic thread adjustment? " 
       << (omp_get_dynamic() ? "yes" : "no")
       << "\nIn parallel? "
       << (omp_in_parallel() ? "yes" : "no")
       << "\nNested parallism? "
       << (omp_get_nested() ? "yes" : "no")
       << "\n" << endl;
#endif
	


  // Parse the input script (each process will just load it's own file
  // for now. ) 

  // Subdomain the grid among the available processors and have each
  // processor set up its grid.

  try {
#ifdef USE_PY_BINDINGS
    if (interactive) {
      PyInterpreter interp(MPI_RANK, MPI_SIZE);
      interp.run();
    } 
#endif

    if (test_run) {
#ifndef HAVE_LIBPOPT
      // Find the exec name
      int argi = 0;
      for (argi = 0; argi < argc; argi++)
      {
        string phred_name(argv[argi]);
        if (phred_name.length() > 5)
          phred_name = phred_name.substr(phred_name.length() - 5, 5);
        
        if (phred_name.compare("phred") == 0)
          break;
      }

      if (argc > (argi + 1))
      {
        string cmd(argv[argi + 1]);

        if (cmd.compare("H") == 0)
          hole();
        else if (cmd.compare("V") == 0 && argi + 4 < argc)
        {
          unsigned int x_cells = 0;
          unsigned int y_cells = 0;
          unsigned int z_cells = 0;

          try
          {
            x_cells = lexical_cast<unsigned int>(argv[argc - 3]);
            y_cells = lexical_cast<unsigned int>(argv[argc - 2]);
            z_cells = lexical_cast<unsigned int>(argv[argc - 1]);
            
            var_benchmark(x_cells, y_cells, z_cells);
          }
          catch(bad_lexical_cast &)
          {
            cout << "To use the variable size benchmark, the last 3 command "
              "line arguments must\nbe natural numbers each greater than "
              "100.\nRunning million node benchmark instead.\n";
            mn_benchmark();
          }
        }
        else if (cmd.compare("M") == 0)
        {
          mn_benchmark();
        }
        else if (cmd.compare("S") == 0 && argi+2 < argc)
        {
          int ysize = 105;
          try
          {
            ysize = lexical_cast<int>(argv[argc - 1]);

            square_hole(ysize);
          }
          catch(bad_lexical_cast &)
          {
            cout << "The size of the hole in the y dimension must "
              "be an integer number of nanometers.\nRunning the "
              "sim with y = 105 nm." << endl;
            square_hole(105);
          }
        } 
        else 
        {
          cout << "Unrecognized test option." << endl;
          usage(0);
        }
        
      } else {
        usage(0);
      }
#else
      mn_benchmark();
#endif
    } else if (!interactive) {

      if (argc > 1)
      {
        string ext = get_extension(argv[argc - 1]);

        // If the file extension is .py, run it in the python interpreter. 
        // If the extension is .jan and load it using Jan's grammer. 

        if (ext.compare("jan") == 0)
        {
          //JanFDTD jfdtd;
          //jfdtd.parse_file(argv[argc - 1]);
          //jfdtd.run(MPI_RANK, MPI_SIZE);
          cout << "Jan's file format is no longer supported." << endl;
        }
        else if (ext.compare("py") == 0)
        {
#ifdef USE_PY_BINDINGS
          PyInterpreter interp(MPI_RANK, MPI_SIZE);
          interp.run_script(argv[argc - 1]);
#else
          cout << "Python support is not compiled into this version." << endl;
#endif
        } else {
          cout << "Unknown input file given.\n\n";
          usage(0);
	}

      } else {

        // TESTS, TEMPORARY
        //point_test(rank, size);
        //pml_test(MPI_RANK, MPI_SIZE);
        //coupler_test(rank,size);
        //takakura_test(rank, size);

        cout << "No filename given to load problem set up from. \n\n";
        usage(0);
      }
    }
  } catch (const std::exception &e) {
    cout << "Caught exception: " << e.what() << endl;
    cout << "Phred terminated with an error. " << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  now = time(NULL);

  int secs = static_cast<int>(now - start);
  int mins = secs / 60;
  secs = secs % 60;
  int hours = mins / 60;
  mins = mins % 60;
  int days = hours / 24;
  hours = hours % 24;

  cout << "Phred is phinished. \nPhred executed for " ;
  
  if (days == 1)
    cout << days << " day, "; 
  else if (days > 1)
    cout << days << " days, "; 
  
  if (hours == 1)
    cout << hours << " hour, ";
  else if (hours > 1)
    cout << hours << " hours, ";
  
  if (mins == 1)
    cout << mins << " minute, ";
  else if (mins > 1)
    cout << mins << " minutes, ";
  
  cout << secs << " seconds." << endl;
  
#ifdef USE_RUSAGE
  struct rusage ru; 
  
  int res = getrusage(RUSAGE_SELF, &ru);
  
  if (res == -1)
  {
    cout << "\nResource usage information is unavailable." << endl;
  }
  else
  {
    cout << "\nResource usage information:";
    cout << "\nMax resident set size: " << ru.ru_maxrss;
    cout << "\nShared text memory size: " << ru.ru_ixrss;
    cout << "\nUnshared data size: " << ru.ru_idrss;
    cout << "\nUnshared stack size: " << ru.ru_isrss;
    cout << "\nPage reclaims: " << ru.ru_minflt;
    cout << "\nPage faults: " << ru.ru_majflt;
    cout << "\nSwaps: " << ru.ru_nswap;
    cout << "\nBlock input operations: " << ru.ru_inblock;
    cout << "\nBlock output operations: " << ru.ru_oublock;
    cout << "\nMessages sent: " << ru.ru_msgsnd;
    cout << "\nMessages recieved: " << ru.ru_msgrcv;
    cout << "\nSignals recieved: " << ru.ru_nsignals;
    cout << "\nVoluntary context switches: " << ru.ru_nvcsw;
    cout << "\nInvoluntary context switches: " << ru.ru_nivcsw;
    cout << endl;
  }

#endif

  // Thank you and goodnight
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  exit (0);
}
