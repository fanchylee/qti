AC_INIT([qti], [0.01], [bug-report@])
AM_INIT_AUTOMAKE([foreign -Wall])
echo "                             Testing for a C compiler"
AC_PROG_CC
echo "                             Testing for a FORTRAN compiler"
AC_PROG_F77
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])

AC_CHECK_HEADERS([gsl/gsl_linalg.h gsl_matrix.h gsl/gsl_vector.h])
AC_HAVE_LIBRARY([m])
AC_HAVE_LIBRARY([gslcblas])
AC_HAVE_LIBRARY([gsl])

AC_OUTPUT
