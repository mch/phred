/* Scanner for Jan's file formats */

/* %option c++ */

%{
#include <stdlib.h>
#include <limits.h>

#include "jan-grammer.hh"
%}

 /* Definitions */

COMMENT ^\*.*
DIGIT [0-9]
NEWLINE \n

%%
 /* Rules */ 

{COMMENT}           yylval.str = yytext; return COMMENT;
{NEWLINE}           return NEWLINE;
[\-+]?[0-9]+("."[0-9]+)?([eE]?[\-+]?[0-9]+)? yylval.dval = strtod(yytext, 0); return NUM;
{DIGIT}+            yylval.ival = strtol(yytext, 0, 10); return INT;
program_mode        yylval.str = yytext; return PROGRAM_MODE;
structure_mode      yylval.str = yytext; return STRUCTURE_MODE;
timestep_mode       yylval.str = yytext; return TIMESTEP_MODE;
runtime             yylval.str = yytext; return RUNTIME;
time_modulo         yylval.str = yytext; return TIME_MODULO;
dimx                yylval.str = yytext; return DIMX;
dimy                yylval.str = yytext; return DIMY;
dimz                yylval.str = yytext; return DIMZ;
deltax              yylval.str = yytext; return DELTAX;
deltay              yylval.str = yytext; return DELTAY;
deltaz              yylval.str = yytext; return DELTAZ;
deltat              yylval.str = yytext; return DELTAT;
nr_of_materials     yylval.str = yytext; return NR_OF_MATERIALS;
[A-Za-z]+           yylval.str = yytext; return STR;

%%
 /* User code */
