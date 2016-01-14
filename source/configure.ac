#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

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

#
#	configure.ac -- input file for configure script
#



AC_PREREQ([2.63])
AC_INIT([FLD-PIPELINE], [0.1], [vigeesh@kis.uni-freiburg.de])
AC_CONFIG_SRCDIR([src])

# Checks for programs.
AC_PROG_CC
m4_include(m4/mpi.m4)
AX_MPI
AC_PROG_CC_C99

set $CC
AC_SUBST(C99_STD,$2)

AC_CHECK_LIB([m], [pow])
AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_LIB([z], [gzread])
AC_CHECK_LIB([curl], [curl_version])

m4_include(m4/cfitsio.m4)
AX_LIB_CFITSIO

m4_include(m4/postgresql.m4)
AX_LIB_POSTGRESQL

m4_include(m4/hdf5.m4)
AX_LIB_HDF5([parallel])

m4_include(m4/drms.m4)
AX_LIB_DRMS

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_MODE_T
AC_C_RESTRICT
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([bzero getcwd memset mkdir pow select sqrt strchr strdup strerror strstr strtol uname])

AC_CONFIG_FILES([Makefile
	src/Makefile])
AC_OUTPUT
echo '-------------------------------------------------------------------------------------'
echo 'Configure: Done'

