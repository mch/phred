/* 
   Phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004-2005 Matt Hughes <mhughe@uvic.ca>

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

#include "PyInterpreter.hh"
#include "../Exceptions.hh"
#include "../Globals.hh"

#include <iostream>
#include <mpi.h>

using namespace std;

#include <boost/python.hpp>
using namespace boost::python;

// Reduce compile time...
void export_results();
void export_excitations();
void export_fdtd();
void export_boundaries();
void export_materials();
void export_types();
void export_datawriters();
void export_grids();
void export_csg();

BOOST_PYTHON_MODULE(Phred)
{
  export_grids();
  export_results();
  export_excitations();
  export_fdtd();
  export_boundaries();
  export_materials();
  export_types();
  export_datawriters();
  export_csg();
}

static struct _inittab modules_[] = 
  {
    {"Phred", &initPhred},
    {0, 0}
  };
  
#include <stdio.h>

#ifdef HAVE_READLINE
#include <unistd.h>
#include <readline/readline.h>
#include <readline/history.h>
#endif

#include <string>

using namespace std;

PyInterpreter::PyInterpreter(int argc, char **argv)
{
  int ret = PyImport_ExtendInittab(modules_);
  if (ret == -1)
    throw PyInterpException("Unable to add Phred modules to Python interpreter.");
  
  Py_Initialize();

  // Setup command line arguments
  setup_args(argc, argv);

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

  add_modules();

  object main_module = extract<object>( PyImport_AddModule("__main__") );
  object main_namespace = main_module.attr("__dict__");

#if PYTHON_VERSION_MAJOR >= 2 && PYTHON_VERSION_MINOR < 3
  PyObject *result = PyRun_File(fp, const_cast<char *>(filename), 
                                Py_file_input, 
                                main_namespace.ptr(), 
                                main_namespace.ptr());
#else
  PyObject *result = PyRun_File(fp, filename, Py_file_input, 
                                main_namespace.ptr(), 
                                main_namespace.ptr());
#endif

  if (!result) {
    PyObject *err = PyErr_Occurred();
    if (err)
    {
      PyErr_Print();
    }
    throw PyInterpException("Python detected an error.");
    //MPI_Abort(MPI_COMM_PHRED, 1);
  }

  fclose(fp);
}

void PyInterpreter::run()
{
  if (MPI_RANK == 0)
    master();
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
  
  add_modules();

  while (size > 0) 
  {
    MPI_Bcast(static_cast<void *>(&size), 1, MPI_INT, 
              0, MPI_COMM_PHRED);

    if (size == -1)
      throw PyInterpException("Standard input and standard output must be bound to a terminal to use interactive mode.");

    // A wasty memory management policy. But it's simple and should
    // always work, barring lack of memory.
    buffer = new char[size + 1];

    if (!buffer)
      throw MemoryException();

    MPI_Bcast(static_cast<void *>(buffer), size, MPI_CHAR, 
              0, MPI_COMM_PHRED);
    buffer[size] = 0;

    try {
      handle<> res( PyRun_String(buffer, Py_file_input, 
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
    MPI_Bcast(static_cast<void *>(&size), 1, MPI_INT, 0, MPI_COMM_PHRED);    

    throw PyInterpException("Standard input and standard output must be bound to a terminal to use interactive mode.");
  }

  handle<> main_module(borrowed( PyImport_AddModule("__main__") ));
  handle<> main_namespace(borrowed( PyModule_GetDict(main_module.get()) ));

  add_modules();

  //if (size_ > 1) {
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
        MPI_Bcast(static_cast<void *>(&size), 1, MPI_INT, 0, MPI_COMM_PHRED);
        MPI_Bcast(static_cast<void *>(const_cast<char *>(buffer.c_str())), buffer.length(), 
                  MPI_CHAR, 0, MPI_COMM_PHRED);

        try {
          handle<> res( PyRun_String(buffer.c_str(), Py_single_input, 
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
    throw PyInterpException("Only one process can be used in interactive mode, and readline is required.");
#endif
    //}  
  // else
//   {
//     char **argv;
//     argv = new char*[2];
//     argv[0] = "phred";
//     argv[1] = "\0";

//     Py_Main(1, argv);  

//     delete[] argv;
//   }
}

void PyInterpreter::setup_args(int argc, char **argv)
{
  object main_module = extract<object>( PyImport_AddModule("__main__") );
  object main_namespace = main_module.attr("__dict__");

  str sys_name("sys");
  object sys_module = extract<object>( PyImport_Import(sys_name.ptr()) );
  dict sys_namespace = extract<dict>( sys_module.attr("__dict__") );

  list args;
  str temp("");
  for (int i = 0; i < argc; i++)
  {
    temp = argv[i];
    args.append( temp );
  }

  sys_namespace["argv"] = args;

  list path = extract<list>( sys_namespace["path"] );
  str dot(".");
  path.insert(0, dot);
}

void PyInterpreter::add_modules()
{
  object main_module = extract<object>( PyImport_AddModule("__main__") );
  dict main_namespace = extract<dict>( main_module.attr("__dict__") );


  // MPI Data
  main_namespace["MPI_RANK"] = MPI_RANK;
  main_namespace["MPI_SIZE"] = MPI_SIZE;

  // Import the contents of the Phred module
  object phred_mod = extract<object>( PyImport_ImportModule("Phred") );
  dict phred_namespace = extract<dict>( phred_mod.attr("__dict__") );
  phred_namespace["__name__"] = "__main__";
  main_namespace.update(phred_namespace);

  // Warrenty text
  str warranty("			    NO WARRANTY\n\
\n\
  11. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY\n\
FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN\n\
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES\n\
PROVIDE THE PROGRAM \"AS IS\" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED\n\
OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF\n\
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS\n\
TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE\n\
PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,\n\
REPAIR OR CORRECTION.\n\
\n\
  12. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING\n\
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR\n\
REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,\n\
INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING\n\
OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED\n\
TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY\n\
YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER\n\
PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE\n\
POSSIBILITY OF SUCH DAMAGES.\n\
");

  main_namespace["warranty"] = warranty;
                       
  // Conditions text
  str conditions("		    GNU GENERAL PUBLIC LICENSE\n\
   TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION\n\
\n\
  0. This License applies to any program or other work which contains\n\
a notice placed by the copyright holder saying it may be distributed\n\
under the terms of this General Public License.  The \"Program\", below,\n\
refers to any such program or work, and a \"work based on the Program\"\n\
means either the Program or any derivative work under copyright law:\n\
that is to say, a work containing the Program or a portion of it,\n\
either verbatim or with modifications and/or translated into another\n\
language.  (Hereinafter, translation is included without limitation in\n\
the term \"modification\".)  Each licensee is addressed as \"you\".\n\
\n\
Activities other than copying, distribution and modification are not\n\
covered by this License; they are outside its scope.  The act of\n\
running the Program is not restricted, and the output from the Program\n\
is covered only if its contents constitute a work based on the\n\
Program (independent of having been made by running the Program).\n\
Whether that is true depends on what the Program does.\n\
\n\
  1. You may copy and distribute verbatim copies of the Program's\n\
source code as you receive it, in any medium, provided that you\n\
conspicuously and appropriately publish on each copy an appropriate\n\
copyright notice and disclaimer of warranty; keep intact all the\n\
notices that refer to this License and to the absence of any warranty;\n\
and give any other recipients of the Program a copy of this License\n\
along with the Program.\n\
\n\
You may charge a fee for the physical act of transferring a copy, and\n\
you may at your option offer warranty protection in exchange for a fee.\n\
\n\
  2. You may modify your copy or copies of the Program or any portion\n\
of it, thus forming a work based on the Program, and copy and\n\
distribute such modifications or work under the terms of Section 1\n\
above, provided that you also meet all of these conditions:\n\
\n\
    a) You must cause the modified files to carry prominent notices\n\
    stating that you changed the files and the date of any change.\n\
\n\
    b) You must cause any work that you distribute or publish, that in\n\
    whole or in part contains or is derived from the Program or any\n\
    part thereof, to be licensed as a whole at no charge to all third\n\
    parties under the terms of this License.\n\
\n\
    c) If the modified program normally reads commands interactively\n\
    when run, you must cause it, when started running for such\n\
    interactive use in the most ordinary way, to print or display an\n\
    announcement including an appropriate copyright notice and a\n\
    notice that there is no warranty (or else, saying that you provide\n\
    a warranty) and that users may redistribute the program under\n\
    these conditions, and telling the user how to view a copy of this\n\
    License.  (Exception: if the Program itself is interactive but\n\
    does not normally print such an announcement, your work based on\n\
    the Program is not required to print an announcement.)\n\
\n\
These requirements apply to the modified work as a whole.  If\n\
identifiable sections of that work are not derived from the Program,\n\
and can be reasonably considered independent and separate works in\n\
themselves, then this License, and its terms, do not apply to those\n\
sections when you distribute them as separate works.  But when you\n\
distribute the same sections as part of a whole which is a work based\n\
on the Program, the distribution of the whole must be on the terms of\n\
this License, whose permissions for other licensees extend to the\n\
entire whole, and thus to each and every part regardless of who wrote it.\n\
\n\
Thus, it is not the intent of this section to claim rights or contest\n\
your rights to work written entirely by you; rather, the intent is to\n\
exercise the right to control the distribution of derivative or\n\
collective works based on the Program.\n\
\n\
In addition, mere aggregation of another work not based on the Program\n\
with the Program (or with a work based on the Program) on a volume of\n\
a storage or distribution medium does not bring the other work under\n\
the scope of this License.\n\
\n\
  3. You may copy and distribute the Program (or a work based on it,\n\
under Section 2) in object code or executable form under the terms of\n\
Sections 1 and 2 above provided that you also do one of the following:\n\
\n\
    a) Accompany it with the complete corresponding machine-readable\n\
    source code, which must be distributed under the terms of Sections\n\
    1 and 2 above on a medium customarily used for software interchange; or,\n\
\n\
    b) Accompany it with a written offer, valid for at least three\n\
    years, to give any third party, for a charge no more than your\n\
    cost of physically performing source distribution, a complete\n\
    machine-readable copy of the corresponding source code, to be\n\
    distributed under the terms of Sections 1 and 2 above on a medium\n\
    customarily used for software interchange; or,\n\
\n\
    c) Accompany it with the information you received as to the offer\n\
    to distribute corresponding source code.  (This alternative is\n\
    allowed only for noncommercial distribution and only if you\n\
    received the program in object code or executable form with such\n\
    an offer, in accord with Subsection b above.)\n\
\n\
The source code for a work means the preferred form of the work for\n\
making modifications to it.  For an executable work, complete source\n\
code means all the source code for all modules it contains, plus any\n\
associated interface definition files, plus the scripts used to\n\
control compilation and installation of the executable.  However, as a\n\
special exception, the source code distributed need not include\n\
anything that is normally distributed (in either source or binary\n\
form) with the major components (compiler, kernel, and so on) of the\n\
operating system on which the executable runs, unless that component\n\
itself accompanies the executable.\n\
\n\
If distribution of executable or object code is made by offering\n\
access to copy from a designated place, then offering equivalent\n\
access to copy the source code from the same place counts as\n\
distribution of the source code, even though third parties are not\n\
compelled to copy the source along with the object code.\n\
\n\
  4. You may not copy, modify, sublicense, or distribute the Program\n\
except as expressly provided under this License.  Any attempt\n\
otherwise to copy, modify, sublicense or distribute the Program is\n\
void, and will automatically terminate your rights under this License.\n\
However, parties who have received copies, or rights, from you under\n\
this License will not have their licenses terminated so long as such\n\
parties remain in full compliance.\n\
\n\
  5. You are not required to accept this License, since you have not\n\
signed it.  However, nothing else grants you permission to modify or\n\
distribute the Program or its derivative works.  These actions are\n\
prohibited by law if you do not accept this License.  Therefore, by\n\
modifying or distributing the Program (or any work based on the\n\
Program), you indicate your acceptance of this License to do so, and\n\
all its terms and conditions for copying, distributing or modifying\n\
the Program or works based on it.\n\
\n\
  6. Each time you redistribute the Program (or any work based on the\n\
Program), the recipient automatically receives a license from the\n\
original licensor to copy, distribute or modify the Program subject to\n\
these terms and conditions.  You may not impose any further\n\
restrictions on the recipients' exercise of the rights granted herein.\n\
You are not responsible for enforcing compliance by third parties to\n\
this License.\n\
\n\
  7. If, as a consequence of a court judgment or allegation of patent\n\
infringement or for any other reason (not limited to patent issues),\n\
conditions are imposed on you (whether by court order, agreement or\n\
otherwise) that contradict the conditions of this License, they do not\n\
excuse you from the conditions of this License.  If you cannot\n\
distribute so as to satisfy simultaneously your obligations under this\n\
License and any other pertinent obligations, then as a consequence you\n\
may not distribute the Program at all.  For example, if a patent\n\
license would not permit royalty-free redistribution of the Program by\n\
all those who receive copies directly or indirectly through you, then\n\
the only way you could satisfy both it and this License would be to\n\
refrain entirely from distribution of the Program.\n\
\n\
If any portion of this section is held invalid or unenforceable under\n\
any particular circumstance, the balance of the section is intended to\n\
apply and the section as a whole is intended to apply in other\n\
circumstances.\n\
\n\
It is not the purpose of this section to induce you to infringe any\n\
patents or other property right claims or to contest validity of any\n\
such claims; this section has the sole purpose of protecting the\n\
integrity of the free software distribution system, which is\n\
implemented by public license practices.  Many people have made\n\
generous contributions to the wide range of software distributed\n\
through that system in reliance on consistent application of that\n\
system; it is up to the author/donor to decide if he or she is willing\n\
to distribute software through any other system and a licensee cannot\n\
impose that choice.\n\
\n\
This section is intended to make thoroughly clear what is believed to\n\
be a consequence of the rest of this License.\n\
\n\
  8. If the distribution and/or use of the Program is restricted in\n\
certain countries either by patents or by copyrighted interfaces, the\n\
original copyright holder who places the Program under this License\n\
may add an explicit geographical distribution limitation excluding\n\
those countries, so that distribution is permitted only in or among\n\
countries not thus excluded.  In such case, this License incorporates\n\
the limitation as if written in the body of this License.\n\
\n\
  9. The Free Software Foundation may publish revised and/or new versions\n\
of the General Public License from time to time.  Such new versions will\n\
be similar in spirit to the present version, but may differ in detail to\n\
address new problems or concerns.\n\
\n\
Each version is given a distinguishing version number.  If the Program\n\
specifies a version number of this License which applies to it and \"any\n\
later version\", you have the option of following the terms and conditions\n\
either of that version or of any later version published by the Free\n\
Software Foundation.  If the Program does not specify a version number of\n\
this License, you may choose any version ever published by the Free Software\n\
Foundation.\n\
\n\
  10. If you wish to incorporate parts of the Program into other free\n\
programs whose distribution conditions are different, write to the author\n\
to ask for permission.  For software which is copyrighted by the Free\n\
Software Foundation, write to the Free Software Foundation; we sometimes\n\
make exceptions for this.  Our decision will be guided by the two goals\n\
of preserving the free status of all derivatives of our free software and\n\
of promoting the sharing and reuse of software generally.\n\
");

  main_namespace["conditions"] = conditions;
}

char *PyInterpreter::rl()
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
