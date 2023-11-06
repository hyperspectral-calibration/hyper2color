#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([hyper2color], [0.5], [ruven@users.sourceforge.net])
AM_INIT_AUTOMAKE


# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor pow])

# Check for OpenMP
AC_OPENMP
AS_IF([test "x$enable_openmp" != "xyes"], [
 INCLUDES="$INCLUDES $OPENMP_CFLAGS"
 LIBS="$LIBS $OPENMP_CFLAGS"
])


AC_CHECK_HEADER([tiff.h],
		[AC_CHECK_LIB(
			[tiff],
			[TIFFOpen],
			[],
			[AC_MSG_ERROR([libtiff not found])]
	)],
	[AC_MSG_ERROR([No libtiff headers found])]
)

# Make sure we can compile and link
AC_MSG_CHECKING([whether libtiff can be compiled])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
			[[#include <tiff.h>]],
			[[int a]])],
		[AC_MSG_RESULT([yes])],
		[AC_MSG_ERROR([unable to compile])]
)

AC_MSG_CHECKING([whether libtiff can be linked])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
			#include <tiffio.h>]],
			[[TIFFGetVersion();]])],
		[AC_MSG_RESULT([yes])],
		[AC_MSG_ERROR([unable to link])]
)


#************************************************************
# Check for Intel Performance Primitives Library

AC_ARG_ENABLE([intel-ipp],
     [  --enable-intel-ipp    enable Intel Integrated Performance Primitives],
     [case "${enableval}" in
       yes) ipp=true ;;
       no)  ipp=false ;;
       *) AC_MSG_ERROR([bad value ${enableval} for --enable-intel-ipp]) ;;
     esac],[ipp=true])
AM_CONDITIONAL([ENABLE_IPP], [test x$ipp = xtrue])

INTEL_IPP=false

if test "x$ipp" = xtrue; then
AC_CHECK_HEADERS( ippi.h,
       AC_SEARCH_LIBS(
               ippiResizeLanczos_8u_C3R,
               ippi,
               INTEL_IPP=true;LIBS="$LIBS -lippcore -lipps -lippcc",
               INTEL_IPP=false,-lippcore -lipps),
	INTEL_IPP=false
)
fi

AM_CONDITIONAL([ENABLE_INTEL_IPP],[test x$INTEL_IPP=xtrue])


AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

AC_CHECK_HEADERS([gsl/gsl_spline.h])

AC_CONFIG_FILES([Makefile \
		 src/Makefile])
AC_OUTPUT
