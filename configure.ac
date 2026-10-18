#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([hyper2color], [0.6], [ruven@users.sourceforge.net])
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
AC_SEARCH_LIBS([floor], [m], [], [AC_MSG_ERROR([libm not found])])

# Check for OpenMP (default on; use --disable-openmp for a serial build)
AC_ARG_ENABLE([openmp],
              [AS_HELP_STRING([--disable-openmp], [build without OpenMP parallelization])],
              [],
              [enable_openmp=yes])
AS_IF([test "x$enable_openmp" = "xyes"], [AC_OPENMP])
AC_SUBST([OPENMP_CFLAGS])


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

AC_SUBST(LIBS)
AC_CONFIG_FILES([Makefile \
		 src/Makefile])
AC_OUTPUT
