# -*- Autoconf -*-


## --------------------------- ##
## Autoconf macros for Numpy. ##
## --------------------------- ##

# CIT_NUMPY_PYTHON_MODULE
# Determine whether the numpy Python module is available.
AC_DEFUN([CIT_NUMPY_PYTHON_MODULE], [
AC_REQUIRE([AM_PATH_PYTHON])
AC_MSG_CHECKING(for numpy python module)
$PYTHON -c "import numpy" 2>/dev/null
if test $? == 0; then
  AC_MSG_RESULT(found)
else
  AC_MSG_FAILURE(not found)
fi
]) dnl CIT_NUMPY_PYTHON_MODULE

# NUMPY_INCDIR
# -----------------
# Determine the directory containing <numpy/arrayobject.h>
AC_DEFUN([CIT_NUMPY_INCDIR], [
AC_REQUIRE([AM_PATH_PYTHON])
AC_CACHE_CHECK([for numpy include directory],
    [_cv_numpy_incdir],
    [_cv_numpy_incdir=`$PYTHON -c "import numpy; numpypath=numpy.__path__[[0]]; print '%s/core/include' % numpypath"`])
AC_SUBST([NUMPY_INCDIR], [$_cv_numpy_incdir])
])dnl CIT_NUMPY_INCDIR


dnl end of file
