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
#   AX_LIB_CFITSIO([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro provides tests of the availability of CFITSIO library.
#
#   The macro adds a --with-cfitsio option accepting the path to the CFITSIO 
#	installation, if the CFITSIO variable is not set. This macro specifically 
#	looks for the fitsio.h file.
#
#   On success, it sets the output variables:
#		CFITSIO_CFLAGS
#		CFITSIO_LDFLAGS
#		CFITSIO_LIBS
#


AU_ALIAS([ACX_LIB_CFITSIO], [AX_LIB_CFITSIO])
AC_DEFUN([AX_LIB_CFITSIO],[

# Default --with-cfitsio configuration option.
AC_ARG_WITH([cfitsio],
  AS_HELP_STRING(
    [--with-cfitsio=[@<:@ARG@:>@]],
            [Path to the CFITSIO Library]
        ),
  [if test "$withval" = "no"; then
     with_cfitsio="no"
   elif test "$withval" = "yes"; then
     with_cfitsio="yes"
   else
     with_cfitsio="yes"
     CFITSIO="$withval"
   fi],
   [with_cfitsio="yes"]
)


CFITSIO_CFLAGS=""
CFITSIO_LDFLAGS=""
CFITSIO_LIBS=""
CFITSIO_VERSION=""

if test "$with_cfitsio" = "yes"; then
	AC_SUBST([CFITSIO])
	if test ! -z "$CFITSIO"; then
		AC_MSG_CHECKING([for fitsio.h])
		AC_COMPILE_IFELSE( [AC_LANG_PROGRAM([[#include <${CFITSIO}/fitsio.h>]],
		              			     [[return 1;]]
		        			   )],
				                   [
				                      AC_MSG_RESULT([yes])
						      		  CFITSIO_CFLAGS=-I$CFITSIO
				                    ],
				                    [
							  			AC_COMPILE_IFELSE( [AC_LANG_PROGRAM([[#include <${CFITSIO}/include/fitsio.h>]],
							  		              			     [[return 1;]]
							  		        			   )],
							  				                   [
							  				                      AC_MSG_RESULT([yes])
							  						      		  CFITSIO_CFLAGS=-I$CFITSIO/include
							  				                    ],
							  				                    [
							  				                      AC_MSG_RESULT([no])
																  echo "Either the CFITSIO path is not set or CFITSIO is NOT installed."
																  AC_MSG_ERROR([Set the CFITSIO path or configure with --with-cfitisio=<PATH>.])
							  				                    ]
							  				                 )
				                    ]
				                 )
		
        CFITSIO_LDFLAGS="-L${CFITSIO}/lib"
        CFITSIO_LIBS="-lcfitsio"
        AC_SUBST(CFITSIO_CFLAGS)
        AC_SUBST(CFITSIO_LDFLAGS)
        AC_SUBST(CFITSIO_LIBS)
	else
	  	echo "Either the CFITSIO path is not set or CFITSIO is NOT installed".
	  	AC_MSG_ERROR([Set the CFITSIO path or configure with --with-cfitisio=<PATH>.])
	fi
fi
])
