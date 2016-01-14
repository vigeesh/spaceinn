# Copyright (C) 2014-2016 Vigeesh Gangadharan <vigeesh@kis.uni-freiburg.de>,
#   Kiepenheuer-Institut fuer Sonnenphysik, Freiburg, Germany.
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

# SYNOPSIS
#
#   AX_LIB_DRMS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro provides tests of the availability of DRMS library.
#
#   The macro adds a --with-drms option accepting the path od the DRMS installation
#	if the DRMS variable is not set.
#
#   On success, it sets the output variables:
#		DRMS_CFLAGS
#		DRMS_LDFLAGS
#		DRMS_LIBS
#

AU_ALIAS([ACX_LIB_DRMS], [AX_LIB_DRMS])
AC_DEFUN([AX_LIB_DRMS],[

# Default --with-drms configuration option.
AC_ARG_WITH([drms],
  AS_HELP_STRING(
    [--with-drms=[@<:@ARG@:>@]],
            [Path to the DRMS Library]
        ),
  [if test "$withval" = "no"; then
     with_drms="no"
   elif test "$withval" = "yes"; then
     with_drms="yes"
   else
     with_drms="yes"
     DRMS="$withval"
   fi],
   [with_drms="yes"]
)

DRMS_CFLAGS=""
DRMS_LDFLAGS=""
DRMS_LIBS=""
DRMS_VERSION=""

if test "$with_drms" = "yes"; then
	AC_SUBST([DRMS])
	if test ! -z "$DRMS"; then
		AC_MSG_CHECKING([for jsoc.h])
		AC_COMPILE_IFELSE( [AC_LANG_PROGRAM([[#include <${DRMS}/include/jsoc.h>]],
		              			     [[return 1;]]
		        			   )],
				                   [
				                  		AC_MSG_RESULT([yes])
								        JSOC_MACH=$(${DRMS}/build/jsoc_machine.csh)
								        DRMS_CFLAGS="-I${DRMS}/include"
								        DRMS_LDFLAGS="-L${DRMS}/lib/${JSOC_MACH}"
								        DRMS_LIBS="-ldrms"
				                    ],
				                    [
										echo ${DRMS}/include/jsoc.h
							  			AC_MSG_RESULT([no])
										echo "Either the DRMS path is not set or DRMS is NOT installed."
										AC_MSG_ERROR([Set the DRMS path or configure with --with-drms=<PATH>.])
				                    ]
				                 )
        AC_SUBST(DRMS_CFLAGS)
        AC_SUBST(DRMS_LDFLAGS)
        AC_SUBST(DRMS_LIBS)
	else
	  	AC_MSG_ERROR([Set the DRMS path or configure with --with-cfitisio=<PATH>.])
	fi
fi
])
