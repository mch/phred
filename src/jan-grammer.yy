%{
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

/* A grammer for loading Jan's FDTD input files */

#include <iostream>
using namespace std;

#include "JanFDTD.hh"
#include "Exceptions.hh"

#define YYERROR_VERBOSE = 1

JanFDTD *janfdtd;
extern FILE *yyin; 

int yylex();
void yyerror(const char *s);
void yyset_in(FILE *in_str);
%}

%union {
  char *str; /* for returning strings */
  double dval; /* for retuinging values */
  int ival; /* for returning integers */
}

%token <dval> JAN_NUM  /* Double precision number */
%token <ival> JAN_INT  /* Integer */
%token <str> JAN_STR  /* String */

%type <str> exprs
%type <dval> exprd
%type <ival> expri
%type <ival> exprm

%token JAN_COMMENT
%token JAN_NEWLINE
%token JAN_PROGRAM_MODE
%token JAN_STRUCTURE_MODE
%token JAN_TIMESTEP_MODE
%token JAN_RUNTIME
%token JAN_TIME_MODULO
%token JAN_DIMX
%token JAN_DIMY
%token JAN_DIMZ
%token JAN_DELTAX
%token JAN_DELTAY
%token JAN_DELTAZ
%token JAN_DELTAT
%token JAN_NR_OF_MATERIALS

%token JAN_ABC_XMIN_TYPE
%token JAN_ABC_YMIN_TYPE
%token JAN_ABC_ZMIN_TYPE
%token JAN_ABC_XMAX_TYPE
%token JAN_ABC_YMAX_TYPE
%token JAN_ABC_ZMAX_TYPE

%token JAN_PML
%token JAN_UPML
%token JAN_EWALL
%token JAN_MWALL
%token JAN_MUR2

%token JAN_SOURCE
%token JAN_WAVEGUIDE
%token JAN_DIPOLE

%token JAN_TYPE
%token JAN_EX
%token JAN_EY
%token JAN_EZ
%token JAN_HX
%token JAN_HY
%token JAN_HZ

%token JAN_FUNCTION
%token JAN_SINE
%token JAN_EXPSINE
%token JAN_GAUSS
%token JAN_GAUSSM 
%token JAN_PULSE

%token JAN_ALPHA
%token JAN_TAU
%token JAN_PERIOD_SIZE
%token JAN_FIXED




%%

input: /* empty */
       | input line
       ;

line: JAN_NEWLINE
      | exprs JAN_NEWLINE { cout << "String Line: " << $1 << endl; }
      | exprd JAN_NEWLINE { cout << "Double Line: " << $1 << endl; }
      | expri JAN_NEWLINE { cout << "Int Line: " << $1 << endl; }
      | exprm JAN_NEWLINE { cout << "Material Line: " << $1 << endl; }
      ;

exprm: expri exprd exprd exprd exprd {
       $$ = $1; janfdtd->add_material($1, $2, $3, $4, $5); }
       ;

exprs: JAN_STR { $$ = $1; }
       | JAN_COMMENT { $$ = ""; }
       | JAN_PROGRAM_MODE JAN_STR { $$ = $2; janfdtd->set_program_mode($2); }
       | JAN_STRUCTURE_MODE JAN_STR { $$ = $2; janfdtd->set_structure_mode($2); }
       | JAN_TIMESTEP_MODE JAN_STR { $$ = $2; janfdtd->set_timestep_mode($2); }
       ;

exprd: JAN_NUM
       | JAN_DELTAX exprd { $$ = $2; janfdtd->set_deltax($2); }
       | JAN_DELTAY exprd { $$ = $2; janfdtd->set_deltay($2); }
       | JAN_DELTAZ exprd { $$ = $2; janfdtd->set_deltaz($2); }
       | JAN_DELTAT exprd { $$ = $2; janfdtd->set_deltat($2); }
       | JAN_DIMX exprd { $$ = $2; janfdtd->set_dimx($2); }
       | JAN_DIMY exprd { $$ = $2; janfdtd->set_dimy($2); }
       | JAN_DIMZ exprd { $$ = $2; janfdtd->set_dimz($2); }
       | JAN_TIME_MODULO expri { $$ = $2; janfdtd->set_time_modulo($2); }
       | JAN_RUNTIME expri { $$ = $2; janfdtd->set_runtime($2); }
       | JAN_NR_OF_MATERIALS exprd { $$ = $2; janfdtd->set_num_materials($2); }
       ;

expri: JAN_INT { $$ = $1; }
       ;

%%

void parse_jan_grammer(const char *filename, JanFDTD *jfdtd)
{
  //yydebug=1;

  janfdtd = jfdtd;
  FILE *fp = fopen(filename, "r");
  /*yyset_in(fp);*/
 yyin = fp;
  yyparse();
}

void yyerror(const char *s)
{
  throw ParserException(s);
}
