# -*- Autoconf -*-


## -------------------------- ##
## Autoconf macros for PETSc. ##
## -------------------------- ##


# CIT_PATH_PETSC([VERSION], [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -----------------------------------------------------------------------
# Check for the PETSc package.  Requires Python.
AC_DEFUN([CIT_PATH_PETSC], [
# $Id: cit_petsc.m4 17942 2011-02-22 20:47:56Z brad $

AC_REQUIRE([AM_PATH_PYTHON])
AC_ARG_VAR(PETSC_DIR, [location of PETSc installation])
AC_ARG_VAR(PETSC_ARCH, [PETSc configuration])

AC_MSG_CHECKING([for PETSc dir])
if test -z "$PETSC_DIR"; then
    AC_MSG_RESULT(no)
    m4_default([$3], [AC_MSG_ERROR([PETSc not found; set PETSC_DIR])])
elif test ! -d "$PETSC_DIR"; then
    AC_MSG_RESULT(no)
    m4_default([$3], [AC_MSG_ERROR([PETSc not found; PETSC_DIR=$PETSC_DIR is invalid])])
elif test ! -d "$PETSC_DIR/include"; then
    AC_MSG_RESULT(broken)
    m4_default([$3], [AC_MSG_ERROR([PETSc include dir $PETSC_DIR/include not found; check PETSC_DIR])])
elif test ! -f "$PETSC_DIR/include/petscversion.h"; then
    AC_MSG_RESULT(broken)
    m4_default([$3], [AC_MSG_ERROR([PETSc header file $PETSC_DIR/include/petscversion.h not found; check PETSC_DIR])])
fi
AC_MSG_RESULT([$PETSC_DIR])

# In what follows, we consistenly check for the new config layout
# first, in case the user is using an old HG working copy with junk in
# it.

AC_MSG_CHECKING([for PETSc arch])
if test -z "$PETSC_ARCH"; then
    if test -d "$PETSC_DIR/conf"; then
        # new config layout; no default config (?)
        AC_MSG_RESULT(no)
        m4_default([$3], [AC_MSG_ERROR([set PETSC_ARCH])])
    elif test ! -f "$PETSC_DIR/bmake/petscconf"; then
        # old config layout (2.3.3 and earlier)
        AC_MSG_RESULT(error)
        m4_default([$3], [AC_MSG_ERROR([PETSc file $PETSC_DIR/bmake/petscconf not found; check PETSC_DIR])])
    else
        cat >petsc.py <<END_OF_PYTHON
[from distutils.sysconfig import parse_makefile

vars = parse_makefile('$PETSC_DIR/bmake/petscconf')
print 'PETSC_ARCH="%s"' % vars['PETSC_ARCH']

]
END_OF_PYTHON
        eval `$PYTHON petsc.py 2>/dev/null`
        rm -f petsc.py
    fi
fi
AC_MSG_RESULT([$PETSC_ARCH])

AC_MSG_CHECKING([for PETSc config])
if test -d "$PETSC_DIR/$PETSC_ARCH/conf"; then
  if test -f "$PETSC_DIR/$PETSC_ARCH/conf/petscvariables"; then
    cit_petsc_petscconf="$PETSC_DIR/$PETSC_ARCH/conf/petscvariables"
  elif test -f "$PETSC_DIR/$PETSC_ARCH/conf/petscconf"; then
    cit_petsc_petscconf="$PETSC_DIR/$PETSC_ARCH/conf/petscconf"
  else 
    AC_MSG_RESULT(no)
    m4_default([$3], [AC_MSG_ERROR([Could not find file with PETSc configuration settings; check PETSC_ARCH/conf])])
  fi
  # Using conf/variables *should* be obsolete for new config.
  #cit_petsc_variables="$PETSC_DIR/conf/variables"
elif test -d "$PESC_DIR/bmake/$PETSC_ARCH"; then
    # old config layout
    cit_petsc_petscconf="$PETSC_DIR/bmake/$PETSC_ARCH/petscconf"
    cit_petsc_variables="$PETSC_DIR/bmake/common/variables"
   if test ! -f "$cit_petsc_variables"; then
       AC_MSG_RESULT(error)
       m4_default([$3], [AC_MSG_ERROR([PETSc config file $cit_petsc_variables not found; check PETSC_DIR])])
   fi
else
    AC_MSG_RESULT(no)
    m4_default([$3], [AC_MSG_ERROR([PETSc config dir not found; check PETSC_ARCH])])
fi
if test ! -f "$cit_petsc_petscconf"; then
    AC_MSG_RESULT(no)
    m4_default([$3], [AC_MSG_ERROR([PETSc config file $cit_petsc_petscconf not found; check PETSC_ARCH])])
fi
AC_MSG_RESULT([$cit_petsc_petscconf])

AC_MSG_CHECKING([for PETSc version == $1])
echo "PETSC_DIR = $PETSC_DIR" > petscconf
echo "PETSC_ARCH = $PETSC_ARCH" >> petscconf
cat $cit_petsc_petscconf $cit_petsc_variables >> petscconf
cat >petsc.py <<END_OF_PYTHON
[from distutils.sysconfig import parse_config_h, parse_makefile, expand_makefile_vars

f = open('$PETSC_DIR/include/petscversion.h')
vars = parse_config_h(f)
f.close()

parse_makefile('petscconf', vars)

keys = (
    'PETSC_VERSION_MAJOR',
    'PETSC_VERSION_MINOR',
    'PETSC_VERSION_SUBMINOR',

    'PETSC_CC_INCLUDES',
    'PETSC_FC_INCLUDES',
    'PETSC_LIB',
    'PETSC_FORTRAN_LIB',

    'CC',
    'CXX',
    'FC',

    'MPI_LIB',
    'MPI_INCLUDE',

    'SIEVE_FLAGS',
)

for key in keys:
    if key[:6] == 'PETSC_':
        value = expand_makefile_vars(str(vars.get(key, '')), vars)
        if key == 'PETSC_LIB':
            # Libtool strips the former.  (Does it ever work?)
            value = value.replace("/System/Library/Frameworks/vecLib.framework/vecLib",
                                  "-Wl,-framework,vecLib")
        print '%s="%s"' % (key, value)
    else:
        print 'PETSC_%s="%s"' % (key, expand_makefile_vars(str(vars.get(key, '')), vars))

]
END_OF_PYTHON
AS_IF([AC_TRY_COMMAND([$PYTHON petsc.py >conftest.sh 2>&AS_MESSAGE_LOG_FD])],
      [],
      [AC_MSG_RESULT(error)
       AC_MSG_FAILURE([cannot parse PETSc configuration])])
eval `cat conftest.sh`
rm -f conftest.sh petsc.py petscconf

[eval `echo $1 | sed 's/\([^.]*\)[.]\([^.]*\)[.]\([^.]*\).*/petsc_1_major=\1; petsc_1_minor=\2; petsc_1_subminor=\3;/'`]
if test -z "$PETSC_VERSION_MAJOR" -o -z "$PETSC_VERSION_MINOR" -o -z "$PETSC_VERSION_SUBMINOR"; then
    AC_MSG_RESULT(no)
    m4_default([$3], [AC_MSG_ERROR([no suitable PETSc package found])])
elif test "$PETSC_VERSION_MAJOR" -eq "$petsc_1_major" -a \
          "$PETSC_VERSION_MINOR" -eq "$petsc_1_minor" -a \
          "$PETSC_VERSION_SUBMINOR" -eq "$petsc_1_subminor" ; then
    AC_MSG_RESULT(yes)
    $2
else
    AC_MSG_RESULT([no ($PETSC_VERSION_MAJOR.$PETSC_VERSION_MINOR.$PETSC_VERSION_SUBMINOR)])
    m4_default([$3], [AC_MSG_ERROR([no suitable PETSc package found])])
fi

AC_SUBST([PETSC_VERSION_MAJOR])
AC_SUBST([PETSC_VERSION_MINOR])
AC_SUBST([PETSC_VERSION_SUBMINOR])
AC_SUBST([PETSC_CC_INCLUDES])
AC_SUBST([PETSC_FC_INCLUDES])
AC_SUBST([PETSC_LIB])
AC_SUBST([PETSC_FORTRAN_LIB])
AC_SUBST([PETSC_CC])
AC_SUBST([PETSC_CXX])
AC_SUBST([PETSC_FC])
AC_SUBST([PETSC_MPI_LIB])
AC_SUBST([PETSC_MPI_INCLUDE])
AC_SUBST([PETSC_SIEVE_FLAGS])
])dnl CIT_PATH_PETSC


# CIT_CHECK_LIB_PETSC
# -------------------
# Try to link against the PETSc libraries.  If the current language is
# C++, determine the value of PETSC_CXX_LIB, which names the extra
# libraries needed when using a C++ compiler.  (As of PETSc v2.3,
# PETSC_CXX_LIB will always be empty; see comment below.)
AC_DEFUN([CIT_CHECK_LIB_PETSC], [
# $Id: cit_petsc.m4 17942 2011-02-22 20:47:56Z brad $
AC_REQUIRE([CIT_PATH_PETSC])dnl
AC_SUBST(PETSC_CXX_LIB)
PETSC_CXX_LIB=
cit_petsc_save_CC=$CC
cit_petsc_save_LIBS=$LIBS
CC=$PETSC_CC
LIBS="$PETSC_LIB $LIBS"
_CIT_LINK_PETSC_IFELSE([], [
    AC_LANG_CASE(
        [C++], [],
        _CIT_CHECK_LIB_PETSC_FAILED
    )
    #
    # Try to guess the correct value for PETSC_CXX_LIB, assuming PETSC_CC
    # is an MPI wrapper.
    #
    # In theory, when PETSC_CC is 'mpicc', *both* the MPI libraries and
    # includes are effectively hidden, and must be extracted in order to
    # use a C++ compiler (the PETSc configuration does not specify a C++
    # compiler command).
    #
    # But this path was only added for symmetry with CIT_HEADER_PETSC.
    # Because, in practice, there is an asymmetry between includes and
    # libs.  When PETSC_CC is 'mpicc', the MPI includes are indeed hidden:
    # PETSC_INCLUDE omits MPI includes.  But PETSC_LIB always explicitly
    # specifies the MPI library, even (redundantly) when PETSC_CC is
    # 'mpicc'.  So, as of PETSc v2.3 at least, this path is never taken.
    CIT_MPI_LIBS(cit_libs, $PETSC_CC, [
	LIBS="$PETSC_LIB $cit_libs $cit_petsc_save_LIBS"
	unset ac_cv_func_PetscInitialize
	_CIT_LINK_PETSC_IFELSE([
	    PETSC_CXX_LIB=$cit_libs
	], [
	    _CIT_CHECK_LIB_PETSC_FAILED
	])
    ], [
	_CIT_CHECK_LIB_PETSC_FAILED
    ])
])
LIBS=$cit_petsc_save_LIBS
CC=$cit_petsc_save_CC
])dnl CIT_CHECK_LIB_PETSC


# _CIT_CHECK_LIB_PETSC_FAILED
# ---------------------------
AC_DEFUN([_CIT_CHECK_LIB_PETSC_FAILED], [
AC_MSG_ERROR([cannot link against PETSc libraries])
])dnl _CIT_CHECK_LIB_PETSC_FAILED


# _CIT_LINK_PETSC_IFELSE([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ----------------------------------------------------------------
AC_DEFUN([_CIT_LINK_PETSC_IFELSE], [
# PetscInitialize() might have C++ linkage.  If the current language
# is C++, allow for this possibility.
AC_LANG_CASE(
    [C++], [
        AC_MSG_CHECKING([for PetscInitialize])
        AC_LINK_IFELSE(_CIT_CHECK_LIB_PETSC_PROGRAM([]), [
                AC_MSG_RESULT([yes (C++)])
                $1
        ], [
            AC_LINK_IFELSE(_CIT_CHECK_LIB_PETSC_PROGRAM([extern "C"]), [
                AC_MSG_RESULT([yes (C)])
                $1
            ], [
                AC_MSG_RESULT(no)
                $2
            ])
        ])
    ],
    [AC_CHECK_FUNC(PetscInitialize, [$1], [$2])]
)
])dnl _CIT_LINK_PETSC_IFELSE


# _CIT_CHECK_LIB_PETSC_PROGRAM
# ----------------------------
AC_DEFUN([_CIT_CHECK_LIB_PETSC_PROGRAM], [
AC_LANG_PROGRAM([[
$1 int PetscInitialize(int *, char ***,const char *,const char *);
]], [[
    PetscInitialize(0, 0, 0, "checklib");
]])
])dnl _CIT_CHECK_LIB_PETSC_PROGRAM


# CIT_HEADER_PETSC
# ----------------
# Try to use PETSc headers.  If the current language is C++, determine
# the value of PETSC_CXX_INCLUDE, which names the extra include paths
# needed when using a C++ compiler... i.e., the MPI includes.  When
# PETSC_CC is set to an MPI wrapper such as 'mpicc', the required MPI
# includes are effectively hidden, and must be extracted in order to
# use a C++ compiler (the PETSc configuration does not specify a C++
# compiler command).
AC_DEFUN([CIT_HEADER_PETSC], [
# $Id: cit_petsc.m4 17942 2011-02-22 20:47:56Z brad $
AC_REQUIRE([CIT_PATH_PETSC])dnl
AC_REQUIRE([CIT_CHECK_LIB_PETSC])dnl
AC_SUBST(PETSC_CXX_INCLUDE)
PETSC_CXX_INCLUDE=
cit_petsc_save_CC=$CC
cit_petsc_save_CPPFLAGS=$CPPFLAGS
cit_petsc_save_LIBS=$LIBS
CC=$PETSC_CC
CPPFLAGS="$PETSC_CC_INCLUDES $CPPFLAGS"
AC_MSG_CHECKING([for petsc.h])
dnl Use AC_TRY_COMPILE instead of AC_CHECK_HEADER because the
dnl latter also preprocesses using $CXXCPP.
AC_TRY_COMPILE([
#include <petsc.h>
], [], [
    AC_MSG_RESULT(yes)
], [
    AC_MSG_RESULT(no)
    AC_LANG_CASE(
        [C++], [],
        _CIT_HEADER_PETSC_FAILED
    )
    # Try to guess the correct value for PETSC_CXX_INCLUDE, assuming
    # PETSC_CC is an MPI wrapper.
    CIT_MPI_INCLUDES(cit_includes, $PETSC_CC, [
	AC_MSG_CHECKING([for petsc.h])
	CPPFLAGS="$PETSC_CC_INCLUDES $cit_includes $cit_petsc_save_CPPFLAGS"
	AC_TRY_COMPILE([
#include <petsc.h>
	], [], [
	    AC_MSG_RESULT(yes)
	    PETSC_CXX_INCLUDE=$cit_includes
	], [
	    AC_MSG_RESULT(no)
	    _CIT_HEADER_PETSC_FAILED
	])
    ], [
	_CIT_HEADER_PETSC_FAILED
    ])
])
AC_LANG_CASE([C++], [
    LIBS="$PETSC_LIB $PETSC_CXX_LIB $LIBS"
    CIT_MPI_CHECK_CXX_LINK(PETSC_CXX_INCLUDE, [$PETSC_LIB],
                           _CIT_TRIVIAL_PETSC_PROGRAM,
                           [whether we can link a trivial C++ PETSc program],
                           [],
			   AC_MSG_FAILURE([cannot link a trivial C++ PETSc program using $CXX]))
])
LIBS=$cit_petsc_save_LIBS
CPPFLAGS=$cit_petsc_save_CPPFLAGS
CC=$cit_petsc_save_CC
])dnl CIT_HEADER_PETSC


# _CIT_HEADER_PETSC_FAILED
# ------------------------
AC_DEFUN([_CIT_HEADER_PETSC_FAILED], [
AC_MSG_ERROR([header "petsc.h" not found])
])dnl _CIT_HEADER_PETSC_FAILED


# _CIT_TRIVIAL_PETSC_PROGRAM
# --------------------------
AC_DEFUN([_CIT_TRIVIAL_PETSC_PROGRAM], [
AC_LANG_PROGRAM([[
#include <petsc.h>
]], [[
    PetscInitialize(0, 0, 0, "trivial");
    PetscFinalize();
]])
])dnl _CIT_TRIVIAL_PETSC_PROGRAM


# CIT_CHECK_LIB_PETSC_SIEVE
# -------------------------
AC_DEFUN([CIT_CHECK_LIB_PETSC_SIEVE], [
AC_MSG_CHECKING([for PETSc/Sieve])
AC_LANG_PUSH(C++)
cit_petsc_save_LIBS=$LIBS
cit_petsc_save_CPPFLAGS=$CPPFLAGS
LIBS="$PETSC_LIB $PETSC_CXX_LIB $LIBS"
CPPFLAGS="$PETSC_CC_INCLUDES $PETSC_CXX_INCLUDE $CPPFLAGS"
AC_LINK_IFELSE(AC_LANG_PROGRAM([[
#include <petscmesh.h>
]], [[
    const int dim = 3;
    ALE::Mesh<int,double> mesh(PETSC_COMM_WORLD, dim);
]]), [
    AC_MSG_RESULT(yes)
], [
    AC_MSG_RESULT(no)
    AC_MSG_FAILURE([cannot build a trivial C++ PETSc program which uses ALE::Sieve])
])
CPPFLAGS=$cit_petsc_save_CPPFLAGS
LIBS=$cit_petsc_save_LIBS
AC_LANG_POP(C++)
])dnl CIT_CHECK_LIB_PETSC_SIEVE


dnl end of file
