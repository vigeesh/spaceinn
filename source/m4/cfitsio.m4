AU_ALIAS([ACX_CFITSIO], [AX_CFITSIO])
AC_DEFUN([AX_CFITSIO],[

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
		AC_COMPILE_IFELSE( [AC_LANG_PROGRAM([[#include <$CFITSIO/fitsio.h>]],
		              			     [[return 1;]]
		        			   )],
				                   [
				                      AC_MSG_RESULT([yes])
						      		  CFITSIO_CFLAGS=-I$CFITSIO
				                    ],
				                    [
							  			AC_COMPILE_IFELSE( [AC_LANG_PROGRAM([[#include <$CFITSIO/include/fitsio.h>]],
							  		              			     [[return 1;]]
							  		        			   )],
							  				                   [
							  				                      AC_MSG_RESULT([yes])
							  						      		  CFITSIO_CFLAGS=-I$CFITSIO/include
							  				                    ],
							  				                    [
							  				                      AC_MSG_RESULT([no])
																  echo Either the CFITSIO path is not set or CFITSIO is NOT installed.
																  AC_MSG_ERROR([Set the CFITSIO path or configure with --with-cfitisio=<PATH>.])
							  				                    ]
							  				                 )
				                    ]
				                 )
		
        #CFITSIO_CFLAGS="-I${CFITSIO}/include"
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
