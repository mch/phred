#include "PyInterpreter.hh"
#include "../Exceptions.hh"

#include <stdio.h>

#ifdef HAVE_READLINE
#include <unistd.h>
#include <readline/readline.h>
#include <readline/history.h>
#endif

#include <string>

using namespace std;

PyInterpreter::PyInterpreter()
{
  int ret = PyImport_ExtendInittab(modules_);
  if (ret == -1)
    throw PyInterpException("Unable to add Phred modules to Python interpreter.");
  
  Py_Initialize();

#ifdef HAVE_READLINE
  rl_bind_key ('\t', rl_insert);
#endif
}

PyInterpreter::~PyInterpreter()
{
  Py_Finalize();
}

void PyInterpreter::run_script(const char *filename)
{
  // Have to use stdio because I need the FILE ptr for Python.
  FILE *fp;
  fp = fopen(filename, "r");

  if (fp == 0)
    throw PyInterpException("Unable to open Python script file.");

  PyRun_SimpleFile(fp, filename);
  fclose(fp);
}

void PyInterpreter::run(int rank, int size)
{
  if (rank == 0)
    maser();
  else
    slave();
}

void PyInterpreter::slave()
{
  // Recieving a message of zero length means that it is time to
  // return. -1 means that the root rank is not bound to a terminal
  // and that we should throw one too. 
  int size = 1;
  char *buffer = 0;

  handle<> main_module(borrowed( PyImport_AddModule("__main__") ));
  handle<> main_namespace(borrowed( PyModule_GetDict(main_module.get()) ));
  
  while (size > 0) 
  {
    MPI_Bcast(static_cast<void *>(&size), 1, MPI_INT, 
              0, MPI_COMM_WORLD);

    if (size == -1)
      throw PyInterpException("Standard input and standard output must be bound to a terminal to use interactive mode.");

    // A wasty memory management policy. But it's simple and should
    // always work, barring lack of memory.
    buffer = new char[size + 1];

    if (!buffer)
      throw MemoryException();

    MPI_Bcast(static_cast<void *>(buffer), size, MPI_CHAR, 
              0, MPI_COMM_WORLD);
    buffer[size] = 0;

    try {
      handle<> res( PyRun_String(buffer.c_str(), Py_file_input, 
                                 main_namespace.get(),
                                 main_namespace.get()));
    }
    catch (error_already_set)
    {
      PyErr_Print();
    }

    delete[] buffer;
  }
}

void PyInterpreter::master()
{
  // If stdin isn't a tty, don't even try. 
  if (!isatty(STDIN_FILENO)) {
    // Tell slaves to abort as well. 
    int size = -1;
    MPI_Bcast(static_cast<void *>(&size), 1, MPI_INT, 0, MPI_COMM_WORLD);    

    throw PyInterpException("Standard input and standard output must be bound to a terminal to use interactive mode.");
  }

  handle<> main_module(borrowed( PyImport_AddModule("__main__") ));
  handle<> main_namespace(borrowed( PyModule_GetDict(main_module.get()) ));
  
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
      buffer += '\n';

      // Ugg, need to check if the last non white space char is a ':', 
      // if so, go multiline, also need to ignore comments. 
      if (strlen(ln) > 0 && ln[strlen(ln) - 1] == ':') 
      {
        prompt = p2;
        multiline = true;
      }

      if ((multiline && strlen(ln) == 0) || (!multiline && buffer.size() > 0))
      {
        int size = buffer.length();
        MPI_Bcast(static_cast<void *>(&size), 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(static_cast<void *>(buffer.c_str()), buffer.length(), 
                  MPI_CHAR, 0, MPI_COMM_WORLD);

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
        buffer.clear();
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

void PyInterpreter::add_modules()
{

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
