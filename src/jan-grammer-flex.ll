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

/* Scanner for Jan's file formats */

/* %option c++ */

%{
#include <stdlib.h>
#include <limits.h>

#include "jan-grammer.h"
%}

 /* Definitions */

COMMENT ^\*.*
DIGIT [0-9]
NEWLINE \n

%%
 /* Rules */ 

{COMMENT}           yylval.str = yytext; return JAN_COMMENT;
{NEWLINE}           return JAN_NEWLINE;
[\-+]?[0-9]+("."[0-9]+)?([eE]?[\-+]?[0-9]+)? yylval.dval = strtod(yytext, 0); return JAN_NUM;
{DIGIT}+            yylval.ival = strtol(yytext, 0, 10); return JAN_INT;
program_mode        yylval.str = yytext; return JAN_PROGRAM_MODE;
structure_mode      yylval.str = yytext; return JAN_STRUCTURE_MODE;
timestep_mode       yylval.str = yytext; return JAN_TIMESTEP_MODE;
runtime             yylval.str = yytext; return JAN_RUNTIME;
time_modulo         yylval.str = yytext; return JAN_TIME_MODULO;
dimx                yylval.str = yytext; return JAN_DIMX;
dimy                yylval.str = yytext; return JAN_DIMY;
dimz                yylval.str = yytext; return JAN_DIMZ;
deltax              yylval.str = yytext; return JAN_DELTAX;
deltay              yylval.str = yytext; return JAN_DELTAY;
deltaz              yylval.str = yytext; return JAN_DELTAZ;
deltat              yylval.str = yytext; return JAN_DELTAT;
nr_of_materials     yylval.str = yytext; return JAN_NR_OF_MATERIALS;
[A-Za-z]+           yylval.str = yytext; return JAN_STR;

%%
 /* User code */
