# -*- Autoconf -*-


## ------------------------ ##
## Autoconf macros for HDF. ##
## ------------------------ ##


# CIT_ARG_HDF5
# ------------
AC_DEFUN([CIT_ARG_HDF5], [
# $Id: cit_hdf.m4 5189 2006-11-07 02:29:23Z leif $
AC_ARG_VAR(PHDF5_HOME, [home path to HDF5 library])
AC_ARG_WITH([hdf5],
    [AC_HELP_STRING([--with-hdf5],
        [enable HDF5 output @<:@default=$1@:>@])],
    [want_hdf5="$withval"],
    [want_hdf5=$1])
])dnl CIT_ARG_HDF5


# CIT_CHECK_LIB_HDF5
# ------------------
AC_DEFUN([CIT_CHECK_LIB_HDF5], [
# $Id: cit_hdf.m4 5189 2006-11-07 02:29:23Z leif $
if test "$want_hdf5" != no; then
    if test -n "$PHDF5_HOME"; then
        LDFLAGS="-L$PHDF5_HOME/lib $LDFLAGS"
    fi
    # check for basic HDF5 function
    AC_SEARCH_LIBS([H5Fopen], [hdf5], [], [
        if test "$want_hdf5" = auto; then
            want_hdf5=no
            AC_MSG_WARN([HDF5 library not found; disabling HDF5 support])
        else
            AC_MSG_ERROR([HDF5 library not found; try setting PHDF5_HOME])
        fi
    ])
fi
])dnl CIT_CHECK_LIB_HDF5


# CIT_CHECK_LIB_HDF5_PARALLEL
# ---------------------------
AC_DEFUN([CIT_CHECK_LIB_HDF5_PARALLEL], [
# $Id: cit_hdf.m4 5189 2006-11-07 02:29:23Z leif $
if test "$want_hdf5" != no; then
    # check for HDF5 parallel-IO function
    AC_SEARCH_LIBS([H5Pset_dxpl_mpio], [hdf5], [], [
        if test "$want_hdf5" = auto; then
            want_hdf5=no
            AC_MSG_WARN([parallel HDF5 library not found; disabling HDF5 support])
        else
            AC_MSG_ERROR([parallel HDF5 library not found; try configuring HDF5 with '--enable-parallel'])
        fi
    ])
fi
])dnl CIT_CHECK_LIB_HDF5_PARALLEL


# CIT_CHECK_HEADER_HDF5
# ---------------------
AC_DEFUN([CIT_CHECK_HEADER_HDF5], [
# $Id: cit_hdf.m4 5189 2006-11-07 02:29:23Z leif $
if test "$want_hdf5" != no; then
    if test -n "$PHDF5_HOME"; then
        CPPFLAGS="-I$PHDF5_HOME/include $CPPFLAGS"
    fi
    AC_CHECK_HEADERS([hdf5.h], [AC_DEFINE([HAVE_HDF5_H])], [
        if test "$want_hdf5" = auto; then
            want_hdf5=no
            AC_MSG_WARN([header 'hdf5.h' not found; disabling HDF5 support])
        else
            AC_MSG_ERROR([header 'hdf5.h' not found])
        fi
    ])
fi
])dnl CIT_CHECK_HEADER_HDF5


dnl end of file
