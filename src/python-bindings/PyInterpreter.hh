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
 */
class PyInterpreter 
{
private:
protected:
public:
  PyInterpreter();
  ~PyInterpreter();

  void run();
};

#endif // PY_INTERPRETER_H
