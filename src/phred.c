/* 
   phred - Phred is a parallel finite difference time domain electromagnetics simulator.

   Copyright (C) 2004 Matt Hughes

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

#include <termios.h>
#include <grp.h>
#include <pwd.h>
*/

#include <stdio.h>
#include <sys/types.h>
#include "system.h"

/* MPI (rocks your socks right off) */
#include <mpi.h>

/* popt plays way nicer with MPI than getopt. Trust me. */
#include <popt.h>

/* Let's use C++ for things that aren't speed critical, because life
   is just so much easier that way. And safer. */
#include <iostream>
#include <vector>
#include <string>

using namespace std; // Too lazy to type namespaces all the time. 

#define EXIT_FAILURE 1

// WTF?
//char *xmalloc ();
//char *xrealloc ();
//char *xstrdup ();


static void usage (int status);

/* The name the program was run with, stripped of any leading path. */
char *program_name;

static struct poptOption options[] = 
  {
    {"help", 'h', POPT_ARG_NON, 0, 'h'},
    {"version", 'V', POPT_ARG_NONE, 0, 'V'},
    {"verbose", 'v', POPT_ARG_NONE, 0, 'v'},
    {"file", 'f', POPT_ARG_STRING, 0 'f'}
  };

static int decode_switches (int argc, char **argv);

int
main (int argc, char **argv)
{
  int i, rank, size;

  MPI_Init(&argc, &args);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0)
    {
      program_name = argv[0];
      
      // MPI implementations are not required to distribute command line
      // args, although MPICH does.
      i = decode_switches (argc, const_cast<const char **>(argv));
    } 
  else 
    {
      program_name = 0;
    }

  /* do the work */
  
  // MPI Goodness?

  // Parse the input script

  // Allocate data structures

  // Main loop
  
  
  // Thank you and goodnight
  MPI_Finalize();

  exit (0);
}

/* Set all the option flags according to the switches specified.
   Return the index of the first non-option argument.  */

static int
decode_switches (int argc, char **argv)
{
  int c;
  char *arg = 0;

  poptContext ctx = poptGetContext(0, argc, argv, options, 0);

  while ((c = poptGetNextOpt (ctx)) > 0 || c == POPT_ERROR_BADOPT)
  {
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

    case 'f':
      arg = const_cast<char *>(poptGetOptArg(ctx));

      if (arg)
        inputfile = arg;
      else 
      {
        cout << "No filename given for the -f switch." << endl;
        exit(0);
      }
      break;

    default:
      cout << "WARNING: got unknown option number: " << c << endl;
    }
  }

  poptFreeContext(ctx);


  
  return optind;
}


static void
usage (int status)
{
  printf (_("%s - \
Phred is a parallel finite difference time domain electromagnetics simulator.\n"), program_name);
  printf (_("Usage: %s [OPTION]... [FILE]...\n"), program_name);
  printf (_("\
Options:\n\
  -i, --interactive          prompt for confirmation\n\
  -v, --verbose              print more information\n\
  -h, --help                 display this help and exit\n\
  -V, --version              output version information and exit\n\
"));
  exit (status);
}
