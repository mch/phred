/* A grammer for loading Jan's FDTD input files */

%{
#include <iostream>
using namespace std;

#define YYERROR_VERBOSE = 1

int yylex();
void yyerror(const char *s);
%}

%union {
  char *str; /* for returning strings */
  double dval; /* for retuinging values */
  int ival; /* for returning integers */
}

%token <dval> NUM  /* Double precision number */
%token <ival> INT  /* Integer */
%token <str> STR  /* String */

%type <str> exprs
%type <dval> exprd
%type <ival> expri
%type <ival> exprm

%token COMMENT
%token NEWLINE
%token PROGRAM_MODE
%token STRUCTURE_MODE
%token TIMESTEP_MODE
%token RUNTIME
%token TIME_MODULO
%token DIMX
%token DIMY
%token DIMZ
%token DELTAX
%token DELTAY
%token DELTAZ
%token DELTAT
%token NR_OF_MATERIALS

%token ABC_XMIN_TYPE
%token ABC_YMIN_TYPE
%token ABC_ZMIN_TYPE
%token ABC_XMAX_TYPE
%token ABC_YMAX_TYPE
%token ABC_ZMAX_TYPE

%token PML
%token EWALL
%token MWALL
%token MUR2

%token SOURCE
%token WAVEGUIDE
%token DIPOLE

%token TYPE
%token EX
%token EY
%token EZ
%token HX
%token HY
%token HZ

%token FUNCTION
%token SINE
%token EXPSINE
%token GAUSS
%token GAUSSM 
%token PULSE

%token ALPHA
%token TAU
%token PERIOD_SIZE
%token FIXED




%%

input: /* empty */
       | input line
       ;

line: NEWLINE
      | exprs NEWLINE { cout << "String Line: " << $1 << endl; }
      | exprd NEWLINE { cout << "Double Line: " << $1 << endl; }
      | expri NEWLINE { cout << "Int Line: " << $1 << endl; }
      | exprm NEWLINE { cout << "Material Line: " << $1 << endl; }
      ;

exprm: exprd exprd exprd exprd exprd {
       $$ = $1; cout << "mat # " << $1 << "eps: "
       << $2 << ", " << $3 << ", etc.\n"; }
       ;

exprs: STR { $$ = $1; }
       | COMMENT { $$ = ""; cout << "parser comment." << endl; }
       | PROGRAM_MODE STR { $$ = $2; cout << "program_mode " <<
       $2 << endl; }
       | STRUCTURE_MODE STR { $$ = $2; cout <<
       "structure_mode " << $2 << endl; }
       | TIMESTEP_MODE STR { $$ = $2; cout <<
       "timestep_mode " << $2 << endl; }
       ;

exprd: NUM
       | DELTAX exprd { $$ = $2; cout << "deltax " <<
       $2 << endl; }
       | DELTAY exprd { $$ = $2; cout << "deltay " <<
       $2 << endl; }
       | DELTAZ exprd { $$ = $2; cout << "deltaz " <<
       $2 << endl; }
       | DELTAT exprd { $$ = $2; cout << "deltat " <<
       $2 << endl; }
       | DIMX exprd { $$ = $2; cout << "dimx " <<
       $2 << endl; }
       | DIMY exprd { $$ = $2; cout << "dimy " <<
       $2 << endl; }
       | DIMZ exprd { $$ = $2; cout << "dimz " <<
       $2 << endl; }
       | TIME_MODULO exprd { $$ = $2; cout << "time_modulo " <<
       $2 << endl; }
       | RUNTIME exprd { $$ = $2; cout << "runtime " <<
       $2 << endl; }
       | NR_OF_MATERIALS exprd { $$ = $2; cout << "num mat: " <<
       $2 << endl; }
       ;

expri: INT { $$ = $1; }
       ;

%%

void parse_jan_grammer()
{
  //yydebug=1;
  yyparse();
}

void yyerror(const char *s)
{
  cout << "Parser error: " << s << endl;
}
