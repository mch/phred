/* 
   phred - Phred is a parallel finite difference time domain
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

#ifndef PY_INTERPRETER_H
#define PY_INTERPRETER_H

#include "../config.h"

/**
 * A wrapper which reads commands from stdin and passes them to a
 * Python interpreter. This is done in this manner because of the need
 * to pass anything the use types to the other processes so they can
 * evaluate them at the same time. All ranks need to stay in lock
 * step. 
 *
 * This functionality will probably only be used for testing, so it
 * would make sense to only allow it when there is only one rank in
 * the communicator, but it's too cool not to do. GNU Readline is
 * used, so if a non GNU licence is needed, we'll restrict it to one
 * rank and use Py_Main().
 *
 * The more important functionality is the loading and execution of
 * user defined script files which set up and solve FDTD problems. 
 */
class PyInterpreter 
{
private:
protected:
  int rank_;
  int size_;

  char *PyInterpreter::rl();
  
public:
  PyInterpreter(int rank, int size);
  ~PyInterpreter();

  /**
   * Runs the interactive interpreter on rank 0 if it is a tty. Throws
   * an exception if it is not. Other ranks run a slave process and
   * MPI is used to keep the slaves in sync with the master.
   */
  void run();

  /**
   * Runs on ranks greater than 0. Recieves commands via MPI from the
   * interactive interpreter running on rank 0 and executes them in it's
   * own python interpreter.
   */
  void slave();

  /**
   * Reads commands from the terminal on rank 0 and broadcasts them to
   * the other ranks before executing the commands itself. 
   */
  void master();

  /**
   * Imports the phred modules into the load namespace so that scripts
   * and interactive users don't have to. Help take some of the
   * mystery out I hope. 
   */
  void add_modules();

  /**
   * Executes a python script file which will normally set up a FDTD
   * problem and solve it. 
   *
   * @param filename the name of the file to execute 
   */
  void run_script(const char *filename);
};

#endif // PY_INTERPRETER_H
