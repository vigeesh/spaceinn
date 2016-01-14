AC_PREREQ([2.63])
AC_INIT([FLD-PIPELINE], [0.1], [vigeesh@kis.uni-freiburg.de])
AC_CONFIG_SRCDIR([src])
#AC_CONFIG_HEADERS([config.h])

#AC_CONFIG_AUX_DIR([bin])

# Checks for programs.
AC_PROG_CC
m4_include(m4/mpi.m4)
AX_MPI

# Checks for libraries.
# Checks for libraries.
# FIXME: Replace `main' with a function in `-lcfitsio':
#AC_CHECK_LIB(cfitsio, ffopen, [cfitsio_present=yes], [cfitsio_present=no], [-L$1])
#if test $cfitsio_present = yes; then
#	#AC_DEFINE(HAVE_CFITSIO,1, [Define if you have the CFITSIO libraries and header files.])
#	CFITSIO_LIBS="-lcfitsio"
#       	CFITSIO_CFLAGS="-I"
#	CFITSIO_LDFLAGS="-L"
#	echo "found CFITSIO"
#else
#	AC_MSG_ERROR([CFITSIO library files not found.])
#fi

m4_include(m4/cfitsio.m4)
AX_CFITSIO

echo $CFITSIO_CFLAGS
echo $CFITSIO_LDFLAGS
echo $CFITSIO_LIBS

# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [pow])
# FIXME: Replace `main' with a function in `-lpq':
#AC_CHECK_LIB([pq], [PQescapeString],
#        [
#        echo "POSTG FOUND"
#        ],
#        [
#        AC_MSG_ERROR([PostGRESQL not found.])
#        ]
#)

m4_include(m4/postgresql.m4)
AX_LIB_POSTGRESQL

echo $POSTGRESQL_CPPFLAGS
echo $POSTGRESQL_LDFLAGS
echo $POSTGRESQL_LIBS
echo $POSTGRESQL_VERSION

# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [pthread_create])
# FIXME: Replace `main' with a function in `-lz':
AC_CHECK_LIB([z], [gzread])

AC_CHECK_LIB(drms,mtrack, [drms_present=yes], [drms_presnt=no])
if test "$drms_present" = "yes"; then
	echo 'DRMS PRESENT.........'
else
	echo 'DRMS NOT PRESENT.........'
fi

#AC_CHECK_LIB(hdf5, H5Fopen, [hdf5_present=yes], [hdf5_present=no], [-L$1])
#    if test $hdf5_present = yes; then
#      echo "hdf found"
#        if test -f /home/vigeesh/hdf5/include/H5config.h; then
#                hdf5par=`grep "HAVE_PARALLEL 1" /home/vigeesh/hdf5/include/H5config.h`
#        elif test -f /home/vigeesh/hdf5/include/H5pubconf.h; then
#                hdf5par=`grep "HAVE_PARALLEL 1" /home/vigeesh/hdf5/include/H5pubconf.h`
#        else
#                echo 'No parallel support'
#        fi
#    else
#        AC_MSG_ERROR([Please install libhdf5 (http://hdf.ncsa.uiuc.edu/HDF5/).])
#    fi

m4_include(m4/hdf5.m4)
AX_LIB_HDF5([parallel])
#if test "$with_hdf5" = "yes"; then
#          MPILD="$HDF5_CC"
#else
#          AC_MSG_ERROR([Unable to find HDF5, we need parallel HDF5.])
#fi

#check for DRMS library
AC_SUBST([DRMS])
if test ! -z "$DRMS"; then
        JSOC_MACH=$(${DRMS}/build/jsoc_machine.csh)
        DRMS_CFLAGS="-I${DRMS}/include"
        DRMS_LDFLAGS="-L${DRMS}/lib/${JSOC_MACH}"
        DRMS_LIBS="-ldrms"
        AC_SUBST(DRMS_CFLAGS)
        AC_SUBST(DRMS_LDFLAGS)
        AC_SUBST(DRMS_LIBS)
else
	AC_MSG_ERROR([Either the DRMS path is not set or NetDRMS is NOT installed.])
fi

echo $DRMS_CFLAGS
echo $DRMS_LDFLAGS
echo $DRMS_LIBS

# Checks for header files.
#AC_CHECK_HEADERS([float.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
#AC_HEADER_STDBOOL
#AC_C_INLINE
#AC_TYPE_MODE_T
#AC_C_RESTRICT
#AC_TYPE_SIZE_T

# Checks for library functions.
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
#AC_CHECK_FUNCS([bzero getcwd memset mkdir pow select sqrt strchr strdup strerror strstr strtol uname])

AC_CONFIG_FILES([Makefile
	src/Makefile])
AC_OUTPUT
echo '-------------------------------------------------------------------------------------'
echo 'Configure: Done'
#echo $CC
#echo $MPICC
#echo $H5CC
#echo $HDF5_CFLAGS

#echo $HDF5_CPPFLAGS
#echo $HDF5_LDFLAGS
#echo $HDF5_LIBS

#echo $CFITSIO_LIBS 
#echo $LIBS
#echo $CFITSIO_CFLAGS

