#include "PyInterpreter.hh"
#include "../Exceptions.hh"

#ifdef HAVE_READLINE
#include <unistd.h>
#include <stdio.h>
#include <readline/readline.h>
#include <readline/history.h>
#endif

#include <string>

using namespace std;

PyInterpreter::PyInterpreter()
{
#ifdef HAVE_READLINE
  rl_bind_key ('\t', rl_insert);
#endif
}

PyInterpreter::~PyInterpreter()
{}

void PyInterpreter::run(int rank, int size)
{
  if (rank == 0)
    maser();
  else
    slave();
}

void PyInterpreter::master()
{
  // If stdin isn't a tty, don't even try. 
  if (!isatty(STDIN_FILENO))
    throw PyInterpException("Standard input and standard output must be bound to a terminal to use interactive mode.");

  handle<> main_module(borrowed( PyImport_AddModule("__main__") ));
  handle<> main_namespace(borrowed( PyModule_GetDict(main_module.get()) ));
  
  handle<> btname ( PyString_FromString("phred") );
  handle<> bt( PyImport_Import(btname.get()) );
  PyDict_SetItemString(main_namespace.get(), "phred", bt.get());
  
  if (size > 1) {
#ifdef HAVE_READLINE
    cout << "Phred interactive Python interpreter running. Type ctrl-d to quit." << endl;
    char *ln = 0; 
    const char *prompt, *p1 = ">>> ", *p2 = "... ";
    prompt = p1;
    bool multiline = false;
    int slen = 0;
    string buffer;
    while ((ln = rl())) {
      buffer += ln;
      
      // Ugg, need to check if the last non white space char is a ':', 
      // if so, go multiline. 
      //if () 
      {
        prompt = p2;
        multiline = true;
      }

      if ((multiline && strlen(ln) == 0) || (!multiline && buffer.size() > 0))
      {
        // MPI_Bcast(buffer.str().c_str(), ...
        try {
          handle<> res( PyRun_String(buffer.c_str(), Py_file_input, 
                                     main_namespace.get(),
                                     main_namespace.get()));
        }
        catch (error_already_set)
        {
          PyErr_Print();
        }
        multiline = false;
        prompt = p1;
      }
    }
    cout << endl;
#else
    throw PyInterpException("Only one process can be used in interactive mode.");
  }  
  else
  {
    char **argv = 0;
    Py_Main(0, argv);  
  }
}

char *PyInterpreter::readline()
{
#ifdef HAVE_READLINE
  char *line_read = readline (">> ");

  /* If the line has any text in it,
     save it on the history. */
  if (line_read && *line_read)
    add_history (line_read);

  return (line_read);

#else
  return 0;
#endif
}
