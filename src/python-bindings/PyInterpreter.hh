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
  struct _inittab modules_[] = 
    {
      {"phred", &initphred},
      {0, 0}
    };
  
public:
  PyInterpreter();
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
