# This file is part of Autoconf.                       -*- Autoconf -*-

# Copyright (C) 2003 Oren Ben-Kiki
# This file is distributed under the same terms as the Autoconf macro files.

# Generate automatic documentation using Doxygen. Works in concert with the
# aminclude.m4 file and a compatible doxygen configuration file. Defines the
# following public macros:
#
# DX_???_FEATURE(ON|OFF) - control the default setting fo a Doxygen feature.
# Supported features are 'DOXYGEN' itself, 'DOT' for generating graphics,
# 'HTML' for plain HTML, 'CHM' for compressed HTML help (for MS users), 'CHI'
# for generating a seperate .chi file by the .chm file, and 'MAN', 'RTF',
# 'XML', 'PDF' and 'PS' for the appropriate output formats. The environment
# variable DOXYGEN_PAPER_SIZE may be specified to override the default 'a4wide'
# paper size.
#
# By default, HTML, PDF and PS documentation is generated as this seems to be
# the most popular and portable combination. MAN pages created by Doxygen are
# usually problematic, though by picking an appropriate subset and doing some
# massaging they might be better than nothing. CHM and RTF are specific for MS
# (note that you can't generate both HTML and CHM at the same time). The XML is
# rather useless unless you apply specialized post-processing to it.
#
# The macro mainly controls the default state of the feature. The use can
# override the default by specifying --enable or --disable. The macros ensure
# that contradictory flags are not given (e.g., --enable-doxygen-html and
# --enable-doxygen-chm, --enable-doxygen-anything with --disable-doxygen, etc.)
# Finally, each feature will be automatically disabled (with a warning) if the
# required programs are missing.
#
# Once all the feature defaults have been specified, call DX_INIT_DOXYGEN with
# the following parameters: a one-word name for the project for use as a
# filename base etc., an optional configuration file name (the default is
# 'Doxyfile', the same as Doxygen's default), and an optional output directory
# name (the default is 'doc').

## ----------##
## Defaults. ##
## ----------##

DX_ENV=""
AC_DEFUN(DX_FEATURE_doc,  ON)
AC_DEFUN(DX_FEATURE_dot,  ON)
AC_DEFUN(DX_FEATURE_man,  OFF)
AC_DEFUN(DX_FEATURE_html, ON)
AC_DEFUN(DX_FEATURE_chm,  OFF)
AC_DEFUN(DX_FEATURE_chi,  OFF)
AC_DEFUN(DX_FEATURE_rtf,  OFF)
AC_DEFUN(DX_FEATURE_xml,  OFF)
AC_DEFUN(DX_FEATURE_pdf,  ON)
AC_DEFUN(DX_FEATURE_ps,   ON)

## --------------- ##
## Private macros. ##
## --------------- ##

# DX_ENV_APPEND(VARIABLE, VALUE)
# ------------------------------
# Append VARIABLE="VALUE" to DX_ENV for invoking doxygen.
AC_DEFUN([DX_ENV_APPEND], [AC_SUBST([DX_ENV], ["$DX_ENV $1='$2'"])])

# DX_DIRNAME_EXPR
# ---------------
# Expand into a shell expression prints the directory part of a path.
AC_DEFUN([DX_DIRNAME_EXPR],
         [[expr ".$1" : '\(\.\)[^/]*$' \| "x$1" : 'x\(.*\)/[^/]*$']])

# DX_IF_FEATURE(FEATURE, IF-ON, IF-OFF)
# -------------------------------------
# Expands according to the M4 (static) status of the feature.
AC_DEFUN([DX_IF_FEATURE], [ifelse(DX_FEATURE_$1, ON, [$2], [$3])])

# DX_REQUIRE_PROG(VARIABLE, PROGRAM)
# ----------------------------------
# Require the specified program to be found for the DX_CURRENT_FEATURE to work.
AC_DEFUN([DX_REQUIRE_PROG], [
AC_PATH_TOOL([$1], [$2])
if test "$DX_FLAG_[]DX_CURRENT_FEATURE$$1" = 1; then
    AC_MSG_WARN([$2 not found - will not DX_CURRENT_DESCRIPTION])
    AC_SUBST([DX_FLAG_[]DX_CURRENT_FEATURE], 0)
fi
])

# DX_TEST_FEATURE(FEATURE)
# ------------------------
# Expand to a shell expression testing whether the feature is active.
AC_DEFUN([DX_TEST_FEATURE], [test "$DX_FLAG_$1" = 1])

# DX_CHECK_DEPEND(REQUIRED_FEATURE, REQUIRED_STATE)
# -------------------------------------------------
# Verify that a required features has the right state before trying to turn on
# the DX_CURRENT_FEATURE.
AC_DEFUN([DX_CHECK_DEPEND], [
test "$DX_FLAG_$1" = "$2" \
|| AC_MSG_ERROR([doxygen-DX_CURRENT_FEATURE ifelse([$2], 1,
                            requires, contradicts) doxygen-DX_CURRENT_FEATURE])
])

# DX_CLEAR_DEPEND(FEATURE, REQUIRED_FEATURE, REQUIRED_STATE)
# ----------------------------------------------------------
# Turn off the DX_CURRENT_FEATURE if the required feature is off.
AC_DEFUN([DX_CLEAR_DEPEND], [
test "$DX_FLAG_$1" = "$2" || AC_SUBST([DX_FLAG_[]DX_CURRENT_FEATURE], 0)
])

# DX_FEATURE_ARG(FEATURE, DESCRIPTION,
#                CHECK_DEPEND, CLEAR_DEPEND,
#                REQUIRE, DO-IF-ON, DO-IF-OFF)
# --------------------------------------------
# Parse the command-line option controlling a feature. CHECK_DEPEND is called
# if the user explicitly turns the feature on (and invokes DX_CHECK_DEPEND),
# otherwise CLEAR_DEPEND is called to turn off the default state if a required
# feature is disabled (using DX_CLEAR_DEPEND). REQUIRE performs additional
# requirement tests (DX_REQUIRE_PROG). Finally, an automake flag is set and
# DO-IF-ON or DO-IF-OFF are called according to the final state of the feature.
AC_DEFUN([DX_ARG_ABLE], [
    AC_DEFUN([DX_CURRENT_FEATURE], [$1])
    AC_DEFUN([DX_CURRENT_DESCRIPTION], [$2])
    AC_ARG_ENABLE(doxygen-$1,
                  [AC_HELP_STRING(DX_IF_FEATURE([$1], [--disable-doxygen-$1],
                                                      [--enable-doxygen-$1]),
                                  DX_IF_FEATURE([$1], [don't $2], [$2]))],
                  [
case "$enableval" in
#(
y|Y|yes|Yes|YES)
    AC_SUBST([DX_FLAG_$1], 1)
    $3
;; #(
n|N|no|No|NO)
    AC_SUBST([DX_FLAG_$1], 0)
;; #(
*)
    AC_MSG_ERROR([invalid value '$enableval' given to doxygen-$1])
;;
esac
], [
AC_SUBST([DX_FLAG_$1], [DX_IF_FEATURE([$1], 1, 0)])
$4
])
if DX_TEST_FEATURE([$1]); then
    $5
    :
fi
if DX_TEST_FEATURE([$1]); then
    AM_CONDITIONAL(DX_COND_$1, :)
    $6
    :
else
    AM_CONDITIONAL(DX_COND_$1, false)
    $7
    :
fi
])

## -------------- ##
## Public macros. ##
## -------------- ##

# DX_XXX_FEATURE(DEFAULT_STATE)
# -----------------------------
AC_DEFUN([DX_DOXYGEN_FEATURE], [AC_DEFUN([DX_FEATURE_doc],  [$1])])
AC_DEFUN([DX_MAN_FEATURE],     [AC_DEFUN([DX_FEATURE_man],  [$1])])
AC_DEFUN([DX_HTML_FEATURE],    [AC_DEFUN([DX_FEATURE_html], [$1])])
AC_DEFUN([DX_CHM_FEATURE],     [AC_DEFUN([DX_FEATURE_chm],  [$1])])
AC_DEFUN([DX_CHI_FEATURE],     [AC_DEFUN([DX_FEATURE_chi],  [$1])])
AC_DEFUN([DX_RTF_FEATURE],     [AC_DEFUN([DX_FEATURE_rtf],  [$1])])
AC_DEFUN([DX_XML_FEATURE],     [AC_DEFUN([DX_FEATURE_xml],  [$1])])
AC_DEFUN([DX_XML_FEATURE],     [AC_DEFUN([DX_FEATURE_xml],  [$1])])
AC_DEFUN([DX_PDF_FEATURE],     [AC_DEFUN([DX_FEATURE_pdf],  [$1])])
AC_DEFUN([DX_PS_FEATURE],      [AC_DEFUN([DX_FEATURE_ps],   [$1])])

# DX_INIT_DOXYGEN(PROJECT, [CONFIG-FILE], [OUTPUT-DOC-DIR])
# ---------------------------------------------------------
# PROJECT also serves as the base name for the documentation files.
# The default CONFIG-FILE is "Doxyfile" and OUTPUT-DOC-DIR is "doc".
AC_DEFUN([DX_INIT_DOXYGEN], [

# Files:
AC_SUBST([DX_PROJECT], [$1])
AC_SUBST([DX_CONFIG], [ifelse([$2], [], Doxyfile, [$2])])
AC_SUBST([DX_DOCDIR], [ifelse([$3], [], doc,      [$3])])

# Environment variables used inside doxygen.cfg:
DX_ENV_APPEND(SRCDIR, $srcdir)
DX_ENV_APPEND(PROJECT, $DX_PROJECT)
DX_ENV_APPEND(DOCDIR, $DX_DOCDIR)
DX_ENV_APPEND(VERSION, $PACKAGE_VERSION)

# Doxygen itself:
DX_ARG_ABLE(doc, [generate documentation],
            [],
            [],
            [DX_REQUIRE_PROG([DX_DOXYGEN], doxygen)
             DX_REQUIRE_PROG([DX_PERL], perl)],
            [DX_ENV_APPEND(PERL_PATH, $DX_PERL)])

# Dot for graphics:
DX_ARG_ABLE(dot, [generate graphics for the documentation],
            [DX_CHECK_DEPEND(doc, 1)],
            [DX_CLEAR_DEPEND(doc, 1)],
            [DX_REQUIRE_PROG([DX_DOT], dot)],
            [DX_ENV_APPEND(HAVE_DOT, YES)
             DX_ENV_APPEND(DOT_PATH, [`DX_DIRNAME_EXPR($DX_DOT)`])],
            [DX_ENV_APPEND(HAVE_DOT, NO)])

# Man pages generation:
DX_ARG_ABLE(man, [generate manual pages],
            [DX_CHECK_DEPEND(doc, 1)],
            [DX_CLEAR_DEPEND(doc, 1)],
            [],
            [DX_ENV_APPEND(GENERATE_MAN, YES)],
            [DX_ENV_APPEND(GENERATE_MAN, NO)])

# RTF file generation:
DX_ARG_ABLE(rtf, [generate RTF documentation],
            [DX_CHECK_DEPEND(doc, 1)],
            [DX_CLEAR_DEPEND(doc, 1)],
            [],
            [DX_ENV_APPEND(GENERATE_RTF, YES)],
            [DX_ENV_APPEND(GENERATE_RTF, NO)])

# XML file generation:
DX_ARG_ABLE(xml, [generate XML documentation],
            [DX_CHECK_DEPEND(doc, 1)],
            [DX_CLEAR_DEPEND(doc, 1)],
            [],
            [DX_ENV_APPEND(GENERATE_XML, YES)],
            [DX_ENV_APPEND(GENERATE_XML, NO)])

# (Compressed) HTML help generation:
DX_ARG_ABLE(chm, [generate compressed HTML help documentation],
            [DX_CHECK_DEPEND(doc, 1)],
            [DX_CLEAR_DEPEND(doc, 1)],
            [DX_REQUIRE_PROG([DX_HHC], hhc)],
            [DX_ENV_APPEND(HHC_PATH, $DX_HHC)
             DX_ENV_APPEND(GENERATE_HTML, YES)
             DX_ENV_APPEND(GENERATE_HTMLHELP, YES)],
            [DX_ENV_APPEND(GENERATE_HTMLHELP, NO)])

# Seperate CHI file generation.
DX_ARG_ABLE(chi, [generate seperate compressed HTML help index file],
            [DX_CHECK_DEPEND(chm, 1)],
            [DX_CLEAR_DEPEND(chm, 1)],
            [],
            [DX_ENV_APPEND(GENERATE_CHI, YES)],
            [DX_ENV_APPEND(GENERATE_CHI, NO)])

# Plain HTML pages generation:
DX_ARG_ABLE(html, [generate plain HTML documentation],
            [DX_CHECK_DEPEND(doc, 1) DX_CHECK_DEPEND(chm, 0)],
            [DX_CLEAR_DEPEND(doc, 1) DX_CLEAR_DEPEND(chm, 0)],
            [],
            [DX_ENV_APPEND(GENERATE_HTML, YES)],
            [DX_TEST_FEATURE(chm) || DX_ENV_APPEND(GENERATE_HTML, NO)])

# PostScript file generation:
DX_ARG_ABLE(ps, [generate PostScript documentation],
            [DX_CHECK_DEPEND(doc, 1)],
            [DX_CLEAR_DEPEND(doc, 1)],
            [DX_REQUIRE_PROG([DX_LATEX], latex)
             DX_REQUIRE_PROG([DX_MAKEINDEX], makeindex)
             DX_REQUIRE_PROG([DX_DVIPS], dvips)
             DX_REQUIRE_PROG([DX_EGREP], egrep)])

# PDF file generation:
DX_ARG_ABLE(pdf, [generate PDF documentation],
            [DX_CHECK_DEPEND(doc, 1)],
            [DX_CLEAR_DEPEND(doc, 1)],
            [DX_REQUIRE_PROG([DX_PDFLATEX], pdflatex)
             DX_REQUIRE_PROG([DX_MAKEINDEX], makeindex)
             DX_REQUIRE_PROG([DX_EGREP], egrep)])

# LaTeX generation for PS and/or PDF:
if DX_TEST_FEATURE(ps) || DX_TEST_FEATURE(pdf); then
    AM_CONDITIONAL(DX_COND_latex, :)
    DX_ENV_APPEND(GENERATE_LATEX, YES)
else
    AM_CONDITIONAL(DX_COND_latex, false)
    DX_ENV_APPEND(GENERATE_LATEX, NO)
fi

# Paper size for PS and/or PDF:
AC_ARG_VAR(DOXYGEN_PAPER_SIZE,
           [a4wide (default), a4, letter, legal or executive])
case "$DOXYGEN_PAPER_SIZE" in
#(
"")
    AC_SUBST(DOXYGEN_PAPER_SIZE, "")
;; #(
a4wide|a4|letter|legal|executive)
    DX_ENV_APPEND(PAPER_SIZE, $DOXYGEN_PAPER_SIZE)
;; #(
*)
    AC_MSG_ERROR([unknown DOXYGEN_PAPER_SIZE='$DOXYGEN_PAPER_SIZE'])
;;
esac

#For debugging:
#echo DX_FLAG_doc=$DX_FLAG_doc
#echo DX_FLAG_dot=$DX_FLAG_dot
#echo DX_FLAG_man=$DX_FLAG_man
#echo DX_FLAG_html=$DX_FLAG_html
#echo DX_FLAG_chm=$DX_FLAG_chm
#echo DX_FLAG_chi=$DX_FLAG_chi
#echo DX_FLAG_rtf=$DX_FLAG_rtf
#echo DX_FLAG_xml=$DX_FLAG_xml
#echo DX_FLAG_pdf=$DX_FLAG_pdf
#echo DX_FLAG_ps=$DX_FLAG_ps
#echo DX_ENV=$DX_ENV
])


# Detect SWIG
# Andrew Collier <colliera@nu.ac.za>
AC_DEFUN([AC_PKG_SWIG],
[
SWIG_REQUEST_VERSION=

changequote(<<, >>)

for a in $1 $2 $3 $4 $5 $6 $7 $8 $9 x; do
    case "$a" in
        x) break;;
        [0-9]*.[0-9]*.[0-9]*) SWIG_REQUEST_VERSION="$a";;
                c++) SWIGFLAGS="$SWIGFLAGS -c++";;
                raw) SWIGFLAGS="$SWIGFLAGS -c";;
    esac
done

changequote([, ])

AC_PATH_PROG(SWIG,swig)

if test -n "$SWIG";
then
        SWIGLIB=`$SWIG -swiglib`

        AC_SUBST(SWIG)
        AC_SUBST(SWIGLIB)
        AC_SUBST(SWIGFLAGS)

        AC_MSG_CHECKING(swig version)

        changequote(<<, >>)
        swig_version=`$SWIG -version 2>&1 | sed 's/.* \([0-9]*\.[0-9]*\.[0-9]*\).*/\1/p; d'`
        swig_major_ver=`expr $swig_version : "\([0-9]\+\)\.[0-9]\+\.*[0-9]*"`
        swig_minor_ver=`expr $swig_version : "[0-9]\+\.\([0-9]\+\)\.*[0-9]*"`
        swig_micro_ver=`expr $swig_version : "[0-9]\+\.[0-9]\+\.*\([0-9]*\)" "|" 0`
        changequote([, ])

        AC_MSG_RESULT($swig_version)

        SWIGVERNUM=`printf "%02d%02d%02d" $swig_major_ver $swig_minor_ver $swig_micro_ver`
        # SWIGVERNUM=`echo $SWIG_REQUEST_VERSION | awk '{ split($[1],a,"\."); print [a[1]*1000000+a[2]*1000+a[3]] }' 2>/dev/null`

        if test -n "$SWIG_REQUEST_VERSION";
        then
                AC_MSG_CHECKING(requested swig version ($SWIG_REQUEST_VERSION))

                changequote(<<, >>)
                swig_major_req=`expr $SWIG_REQUEST_VERSION : '\([0-9]*\)\.[0-9]*\.[0-9]*'`
                swig_minor_req=`expr $SWIG_REQUEST_VERSION : '[0-9]*\.\([0-9]*\)\.[0-9]*'`
                swig_micro_req=`expr $SWIG_REQUEST_VERSION : '[0-9]*\.[0-9]*\.\([0-9]*\)'`
                changequote([, ])

                if test $swig_major_ver -ge $swig_major_req &&
                   test $swig_minor_ver -ge $swig_minor_req &&
                   test $swig_micro_ver -ge $swig_micro_req
                then
                        AC_MSG_RESULT(yes)
                else
                        AC_MSG_RESULT(no)
                        AC_MSG_ERROR(installed version of swig is too old!)
                fi
        fi
fi
])

dnl @synopsis AC_LIB_WAD
dnl
dnl This macro searches for installed WAD library.
dnl
AC_DEFUN([AC_LIB_WAD],
[
AC_ARG_ENABLE(wad,
        AC_HELP_STRING([--enable-wad], [enable wad module]),
        [
                case "${enableval}" in
                no)  ;;
            *)   if test "x${enableval}" = xyes;
                                 then
                                        check_wad="yes"
                                 fi ;;
        esac
        ], [])

if test -n "$check_wad";
then
        # this won't work unless PYTHON_LINK and PYTHON_EXTRA_LIBS defined
        AC_CHECK_LIB(wadpy, _init, [WADPY=-lwadpy], [], $PYTHON_LINK $PYTHON_EXTRA_LIBS)
        AC_SUBST(WADPY)
fi
])

dnl ----------------------------------------------------------------------------
dnl @synopsis AC_CXX_LIB_BLITZ([optional-string "required"])
dnl
dnl Check whether Blitz++ is installed.
dnl Blitz++ is available at http://oonumerics.org/blitz.
dnl
dnl Set the path for Blitz++ with the option
dnl --with-blitz[=DIR]
dnl Blitz headers should be under DIR/includes
dnl Blitz library should be under DIR/lib
dnl Then try to compile and run a simple program with a Blitz Array
dnl Optional argument `required' triggers an error if Blitz++ not installed
dnl dnl @version $Id: ac_cxx_lib_blitz.m4,v 1.4 2004/03/25 14:17:52 patricg Exp $
dnl @author Patrick Guio <address@bogus.example.com>
dnl
AC_DEFUN([AC_MSG_ERROR_BLITZ],[
AC_MSG_ERROR([
$PACKAGE_STRING requires the Blitz++ template library
available at http://oonumerics.org/blitz
When installed give the directory of installation with the option
--with-blitz@<:@=DIR@:>@
])])


AC_DEFUN([AC_CXX_LIB_BLITZ],[


AC_ARG_WITH(blitz,
AC_HELP_STRING([--with-blitz@<:@=DIR@:>@],[Set the path for Blitz++]),
[],[withval='yes'])


if test "$1" = required -a "$withval" = no ; then
        AC_MSG_ERROR_BLITZ
fi


if test "$withval" != no ; then


        saveCPPFLAGS=$CPPFLAGS
        saveLDFLAGS=$LDFLAGS
        saveLIBS=$LIBS


        if test "$withval" != 'yes'; then
                CPPFLAGS="-I$withval/include"
                LDFLAGS="-L$withval/lib"
        fi
        LIBS="-lblitz"


        AC_CACHE_CHECK([whether Blitz++ is installed],ac_cxx_lib_blitz,
        [AC_LANG_SAVE
        AC_LANG_CPLUSPLUS
        AC_RUN_IFELSE(
        [AC_LANG_PROGRAM([[
#include <blitz/array.h>
]],[[
blitz::Array<int,1> x(10);
x = blitz::tensor::i;
        ]])],[ac_cxx_lib_blitz=yes],[ac_cxx_lib_blitz=no])
        AC_LANG_RESTORE
        ])


        CPPFLAGS=$saveCPPFLAGS
        LDFLAGS=$saveLDFLAGS
        LIBS=$saveLIBS


        if test "$ac_cxx_lib_blitz" = yes ; then
                if test "$withval" != yes ; then
                        CPPFLAGS="-I$withval/include $CPPFLAGS"
                        LDFLAGS="-L$withval/lib $LDFLAGS"
                fi
                LIBS="-lblitz $LIBS"
                AC_DEFINE([HAVE_BLITZ], [1], [Using Blitz++])
        else
                if test "$1" = required ; then
                        AC_MSG_ERROR_BLITZ
                fi
        fi


fi


])


dnl -------------------------------------------------------------------------
dnl @synopsis AC_CXX_PYTHON([optional-string "required"])
dnl
dnl Check whether Python is installed.
dnl  Python is available at http://python.org/.
dnl
AC_DEFUN([AC_CXX_PYTHON], [

AC_ARG_WITH(python, AC_HELP_STRING([--with-python], 
                    [Generate Python scripting support]))
AC_ARG_WITH(python_libs, AC_HELP_STRING([--with-python-libs], 
                         [Location of Python libraries]))
AC_ARG_WITH(python_includes, AC_HELP_STRING([--with-python-includes], 
                             [Location of Python include files]))

if test "$with_python" != no ; then

PYTHON_LDFLAGS=`echo "import distutils.sysconfig; print distutils.sysconfig.get_config_var('LINKFORSHARED')" | python -`
PYTHON_INCLUDES=`echo "import distutils.sysconfig; print distutils.sysconfig.get_config_var('INCLUDEPY')" | python -`
PYTHON_VERSION=`echo "import distutils.sysconfig; print distutils.sysconfig.get_config_var('VERSION')" | python -`

        saveCPPFLAGS=$CPPFLAGS
        saveCXXFLAGS=$CXXFLAGS
        saveLDFLAGS=$LDFLAGS
        saveLIBS=$LIBS

        CPPFLAGS="$CPPFLAGS -I$PYTHON_INCLUDES"

        if [[ ! -z "$with_python_includes" ]]; then 
          CPPFLAGS="$CPPFLAGS -I$with_python_includes"
        fi

        AC_CHECK_HEADER([Python.h], [], [])

        if [[ "$target_vendor" = "apple" ]]; then       
                LIBS="$LIBS -framework Python"
        else
                LIBS="$LIBS -lpython$PYTHON_VERSION"
        fi

        AC_CACHE_CHECK([whether Python is installed],ac_cxx_python,
        [AC_LANG_SAVE
        AC_LANG_CPLUSPLUS
        AC_RUN_IFELSE(
        [AC_LANG_PROGRAM([[
#include <Python.h>
]],[[ Py_Initialize();
        ]])],[ac_cxx_python=yes],[ac_cxx_python=no])
        AC_LANG_RESTORE
        ])

        CPPFLAGS=$saveCPPFLAGS
        CXXFLAGS=$saveCXXFLAGS
        LDFLAGS=$saveLDFLAGS
        LIBS=$saveLIBS

        if test "$ac_cxx_python" = "yes" ; then
                AC_DEFINE([HAVE_PYTHON], [1], [Using Python])
                LDFLAGS="$LDFLAGS $PYTHON_LDFLAGS"
                CPPFLAGS="$CPPFLAGS -I$PYTHON_INCLUDES"
                if [[ ! -z "$with_python_includes" ]]; then 
                        CPPFLAGS="$CPPFLAGS -I$with_python_includes"
                fi
                if [[ "$target_vendor" = "apple" ]]; then       
                        LIBS="$LIBS -framework Python"
                else
                        LIBS="$LIBS -lpython$PYTHON_VERSION"
                fi

        fi
        
fi

])

dnl -------------------------------------------------------------------------
dnl @synopsis AC_CXX_LIB_BOOST_PYTHON([optional-string "required"])
dnl
dnl Check whether Boost Python is installed.
dnl Boost Python is available at http://boost.org/.
dnl
dnl Set the path to the Boost installation with
dnl --with-boost[=DIR]
dnl
dnl Then try to compile and run a simple program with a Blitz Array
dnl Optional argument `required' triggers an error if Blitz++ not installed
dnl dnl @version $Id: ac_cxx_lib_boost_python.m4,v 1.4 2004/03/25 14:17:52 patricg Exp $
dnl @author Patrick Guio <address@bogus.example.com>
dnl
AC_DEFUN([AC_MSG_ERROR_BOOST_PYTHON],[
AC_MSG_ERROR([
$PACKAGE_STRING requires the Boost Python library
available at http://boost.org/
Give the location of the boost installation with
--with-boost@<:@=DIR@:>@
])])


AC_DEFUN([AC_CXX_LIB_BOOST_PYTHON],[


AC_ARG_WITH(boost,
AC_HELP_STRING([--with-boost@<:@=DIR@:>@],[Set the path to the Boost installation]),
[],[withval='yes'])

if test "$1" = required -a "$withval" = no ; then
        AC_MSG_ERROR_BOOST_PYTHON
fi

if test "$withval" != no ; then

        saveCPPFLAGS=$CPPFLAGS
        saveLDFLAGS=$LDFLAGS
        saveLIBS=$LIBS

        if [[ "$target_vendor" = "apple" ]]; then       
                LIBS="$LIBS -framework Python"
        fi

        if test "$withval" != 'yes'; then
                CPPFLAGS="-I$withval/include"
                LDFLAGS="-L$withval/lib"
        fi
        LIBS="$LIBS -lboost_python"


        AC_CACHE_CHECK([whether Boost Python is installed],ac_cxx_lib_boost_python,
        [AC_LANG_SAVE
        AC_LANG_CPLUSPLUS
        AC_RUN_IFELSE(
        [AC_LANG_PROGRAM([[
#include <Python.h>
#include <boost/python.hpp>
char const* greet()
{ return "hello, world"; }
using namespace boost::python;
BOOST_PYTHON_MODULE(hello)
{ def("greet", greet); }
]],[[
        ]])],[ac_cxx_lib_boost_python=yes],[ac_cxx_lib_boost_python=no])
        AC_LANG_RESTORE
        ])


        CPPFLAGS=$saveCPPFLAGS
        LDFLAGS=$saveLDFLAGS
        LIBS=$saveLIBS


        if test "$ac_cxx_lib_boost_python" = yes ; then
                if test "$withval" != yes ; then
                        CPPFLAGS="-I$withval/include $CPPFLAGS"
                        LDFLAGS="-L$withval/lib $LDFLAGS"
                fi
                LIBS="-lboost_python $LIBS"
                AC_DEFINE([HAVE_BOOST_PYTHON], [1], [Using Boost Python])
        else
                if test "$1" = required ; then
                        AC_MSG_ERROR_BOOST_PYTHON
                fi
        fi


fi


])



dnl -------------------------------------------------------------------------
dnl @synopsis AC_CXX_LIB_VTK([optional-string "required"])
dnl
dnl Check whether VTK is installed.
dnl VTK is available at http://public.kitware.com/VTK/.
dnl
dnl Set the path for VTK with the option
dnl --with-vtk[=DIR]
dnl VTK headers should be under DIR/includes
dnl VTK library should be under DIR/lib
dnl Then try to compile and run a simple program with a VTK Sphere
dnl Optional argument `required' triggers an error if VTK not installed
dnl @author Matt Hughes <address@bogus.example.com>
dnl
AC_DEFUN([AC_MSG_ERROR_VTK],[
AC_MSG_ERROR([
$PACKAGE_STRING requires the VTK library
available at http://public.kitware.com/VTK/
When installed give the directory of installation with the option
--with-vtk@<:@=DIR@:>@
])])


AC_DEFUN([AC_CXX_LIB_VTK],[


AC_ARG_WITH(vtk,
AC_HELP_STRING([--with-vtk@<:@=DIR@:>@],[Set the path for VTK]),
[],[withval='yes'])


if test "$1" = required -a "$withval" = no ; then
        AC_MSG_ERROR_VTK
fi


if test "$withval" != no ; then


        saveCPPFLAGS=$CPPFLAGS
        saveLDFLAGS=$LDFLAGS
        saveLIBS=$LIBS


        if test "$withval" != 'yes'; then
                CPPFLAGS="-I$withval/include"
                LDFLAGS="-L$withval/lib"
        fi
        LIBS="-lvtk"


        AC_CACHE_CHECK([whether VTK is installed],ac_cxx_lib_vtk,
        [AC_LANG_SAVE
        AC_LANG_CPLUSPLUS
        AC_RUN_IFELSE(
        [AC_LANG_PROGRAM([[
#include "vtkStructuredGrid.h"
]],[[
  static int dims[3]={13,11,11};
  
  // Create the structured grid.
  vtkStructuredGrid *sgrid = vtkStructuredGrid::New();
  sgrid->SetDimensions(dims);

        ]])],[ac_cxx_lib_vtk=yes],[ac_cxx_lib_vtk=no])
        AC_LANG_RESTORE
        ])


        CPPFLAGS=$saveCPPFLAGS
        LDFLAGS=$saveLDFLAGS
        LIBS=$saveLIBS


        if test "$ac_cxx_lib_vtk" = yes ; then
                if test "$withval" != yes ; then
                        CPPFLAGS="-I$withval/include $CPPFLAGS"
                        LDFLAGS="-L$withval/lib $LDFLAGS"
                fi
                LIBS="-lblitz $LIBS"
                AC_DEFINE([HAVE_VTK], [1], [Using VTK])
        else
                if test "$1" = required ; then
                        AC_MSG_ERROR_VTK
                fi
        fi


fi


])

