# -*- Autoconf -*-


## ------------------------ ##
## Autoconf macros for MPI. ##
## ------------------------ ##


# CIT_PROG_MPICC
# --------------
# Call AC_PROG_CC, but prefer MPI C wrappers to a bare compiler in
# the search list.  Set MPICC to the program/wrapper used to compile
# C MPI programs.  Set CC to the compiler used to compile ordinary
# C programs, and link shared libraries of all types (see the
# comment about the MPI library, below).  Make sure that CC and
# MPICC both represent the same underlying C compiler.
AC_DEFUN([CIT_PROG_MPICC], [
# $Id: cit_mpi.m4 7980 2007-09-18 01:01:59Z leif $
AC_PROVIDE([_CIT_PROG_MPICC])dnl
AC_REQUIRE([_CIT_PROG_MPICC_SEARCH_LIST])dnl
AC_BEFORE([$0], [AC_PROG_CC])
AC_ARG_VAR(MPICC, [MPI C compiler command])
AC_SUBST([MPICC])
test -z "$want_mpi" && want_mpi=yes
# The 'cit_compiler_search_list' is the result of merging the
# following:
#     * MPI C wrappers
#     * the range of values for config's COMPILER_CC_NAME
#       (cc cl ecc gcc icc pgcc xlc xlc_r)
# Newer names are tried first (e.g., icc before ecc).
cit_compiler_search_list="gcc cc cl icc ecc pgcc xlc xlc_r"
# There are two C command variables, so there are four cases to
# consider:
#
#     ./configure CC=gcc MPICC=mpicc       # save MPICC as cit_MPICC; MPICC=$CC
#     ./configure CC=gcc                   # MPICC=$CC, guess cit_MPICC
#     ./configure MPICC=mpicc              # derive CC
#     ./configure                          # guess MPICC and derive CC
#
# In the cases where CC is explicitly specified, the MPI C wrapper
# (cit_MPICC, if known) is only used to gather compile/link flags (if
# needed).
if test "$want_mpi" = yes; then
    if test -n "$CC"; then
        cit_MPICC_underlying_CC=$CC
        if test -n "$MPICC"; then
            # CC=gcc MPICC=mpicc
            cit_MPICC=$MPICC
            MPICC=$CC
        else
            # CC=gcc MPICC=???
            AC_CHECK_PROGS(cit_MPICC, $cit_mpicc_search_list)
        fi
    else
        if test -n "$MPICC"; then
            # CC=??? MPICC=mpicc
            cit_MPICC=$MPICC
            CC=$MPICC # will be reevaluated below
        else
            # CC=??? MPICC=???
            cit_compiler_search_list="$cit_mpicc_search_list $cit_compiler_search_list"
        fi
    fi
fi
AC_PROG_CC($cit_compiler_search_list)
if test "$want_mpi" = yes; then
    if test -z "$MPICC"; then
        MPICC=$CC
    fi
    if test -z "$cit_MPICC"; then
        case $MPICC in
            *mp* | hcc)
                cit_MPICC=$MPICC
                ;;
        esac
    fi
    # The MPI library is typically static.  Linking a shared object
    # against static library is non-portable, and needlessly bloats our
    # Python extension modules on the platforms where it does work.
    # Unless CC was set explicitly, attempt to set CC to the underlying
    # compiler command, so that we may link with the matching C
    # compiler, but omit -lmpi/-lmpich from the link line.
    if test -z "$cit_MPICC_underlying_CC"; then
        if test -n "$cit_MPICC"; then
            AC_MSG_CHECKING([for the C compiler underlying $cit_MPICC])
            CC=
            AC_LANG_PUSH(C)
            for cit_arg_show in CIT_MPI_COMPILE_INFO_SWITCHES
            do
                cit_cmd="$cit_MPICC -c $cit_arg_show"
                if $cit_cmd >/dev/null 2>&1; then
                    CC=`$cit_cmd 2>/dev/null | sed 's/ .*//'`
                    if test -n "$CC"; then
                        AC_COMPILE_IFELSE([AC_LANG_PROGRAM()], [break 2], [CC=])
                    fi
                fi
            done
            AC_LANG_POP(C)
            if test -n "$CC"; then
                AC_MSG_RESULT($CC)
            else
                AC_MSG_RESULT(failed)
                AC_MSG_FAILURE([can not determine the C compiler underlying $cit_MPICC])
            fi
        fi
        cit_MPICC_underlying_CC=$CC
    fi
fi
])dnl CIT_PROG_MPICC


# _CIT_PROG_MPICC
# ---------------
# Search for an MPI C wrapper. ~ This private macro is employed by
# C++-only projects (via CIT_CHECK_LIB_MPI and CIT_HEADER_MPI).  It
# handles the case where an MPI C wrapper is present, but an MPI C++
# wrapper is missing or broken.  This can happen if a C++ compiler was
# not found/specified when MPI was installed.
AC_DEFUN([_CIT_PROG_MPICC], [
AC_REQUIRE([_CIT_PROG_MPICC_SEARCH_LIST])dnl
AC_CHECK_PROGS(cit_MPICC, $cit_mpicc_search_list)
])dnl _CIT_PROG_MPICC


# _CIT_PROG_MPICC_SEARCH_LIST
# ---------------------------
AC_DEFUN([_CIT_PROG_MPICC_SEARCH_LIST], [
# $Id: cit_mpi.m4 7980 2007-09-18 01:01:59Z leif $
cit_mpicc_search_list="mpicc hcc mpcc mpcc_r mpxlc cmpicc"
])dnl _CIT_PROG_MPICC_SEARCH_LIST


# CIT_PROG_MPICXX
# ---------------
# Call AC_PROG_CXX, but prefer MPI C++ wrappers to a bare compiler in
# the search list.  Set MPICXX to the program/wrapper used to compile
# C++ MPI programs.  Set CXX to the compiler used to compile ordinary
# C++ programs, and link shared libraries of all types (see the
# comment about the MPI library, below).  Make sure that CXX and
# MPICXX both represent the same underlying C++ compiler.
AC_DEFUN([CIT_PROG_MPICXX], [
# $Id: cit_mpi.m4 7980 2007-09-18 01:01:59Z leif $
AC_BEFORE([$0], [AC_PROG_CXX])
AC_ARG_VAR(MPICXX, [MPI C++ compiler command])
AC_SUBST([MPICXX])
test -z "$want_mpi" && want_mpi=yes
# The 'cit_compiler_search_list' is the result of merging the
# following:
#     * MPI C++ wrappers
#     * the Autoconf default (g++ c++ gpp aCC CC cxx cc++ cl
#       FCC KCC RCC xlC_r xlC)
#     * the range of values for config's COMPILER_CXX_NAME (aCC CC cl
#       cxx ecpc g++ icpc KCC pgCC xlC xlc++_r xlC_r)
# Newer names are tried first (e.g., icpc before ecpc).
cit_compiler_search_list="g++ c++ gpp aCC CC cxx cc++ cl FCC KCC RCC xlc++_r xlC_r xlC"
cit_compiler_search_list="$cit_compiler_search_list icpc ecpc pgCC"
cit_mpicxx_search_list="mpicxx mpic++ mpiCC hcp mpCC mpxlC mpxlC_r cmpic++"
# There are two C++ command variables, so there are four cases to
# consider:
#
#     ./configure CXX=g++ MPICXX=mpicxx    # save MPICXX as cit_MPICXX; MPICXX=$CXX
#     ./configure CXX=g++                  # MPICXX=$CXX, guess cit_MPICXX
#     ./configure MPICXX=mpicxx            # derive CXX
#     ./configure                          # guess MPICXX and derive CXX
#
# In the cases where CXX is explicitly specified, the MPI C++ wrapper
# (cit_MPICXX, if known) is only used to gather compile/link flags (if
# needed).
if test "$want_mpi" = yes; then
    if test -n "$CXX"; then
        cit_MPICXX_underlying_CXX=$CXX
        if test -n "$MPICXX"; then
            # CXX=g++ MPICXX=mpicxx
            cit_MPICXX=$MPICXX
            MPICXX=$CXX
        else
            # CXX=g++ MPICXX=???
            AC_CHECK_PROGS(cit_MPICXX, $cit_mpicxx_search_list)
        fi
    else
        if test -n "$MPICXX"; then
            # CXX=??? MPICXX=mpicxx
            cit_MPICXX=$MPICXX
            CXX=$MPICXX # will be reevaluated below
        else
            # CXX=??? MPICXX=???
            cit_compiler_search_list="$cit_mpicxx_search_list $cit_compiler_search_list"
        fi
    fi
fi
AC_PROG_CXX($cit_compiler_search_list)
if test "$want_mpi" = yes; then
    if test -z "$MPICXX"; then
        MPICXX=$CXX
    fi
    if test -z "$cit_MPICXX"; then
        case $MPICXX in
            *mp* | hcp)
                cit_MPICXX=$MPICXX
                ;;
        esac
    fi
    # The MPI library is typically static.  Linking a shared object
    # against static library is non-portable, and needlessly bloats our
    # Python extension modules on the platforms where it does work.
    # Unless CXX was set explicitly, attempt to set CXX to the underlying
    # compiler command, so that we may link with the matching C++
    # compiler, but omit -lmpi/-lmpich from the link line.
    if test -z "$cit_MPICXX_underlying_CXX"; then
        if test -n "$cit_MPICXX"; then
            AC_MSG_CHECKING([for the C++ compiler underlying $cit_MPICXX])
            CXX=
            AC_LANG_PUSH(C++)
            for cit_arg_show in CIT_MPI_COMPILE_INFO_SWITCHES
            do
                cit_cmd="$cit_MPICXX -c $cit_arg_show"
                if $cit_cmd >/dev/null 2>&1; then
                    CXX=`$cit_cmd 2>/dev/null | sed 's/ .*//'`
                    if test -n "$CXX"; then
                        AC_COMPILE_IFELSE([AC_LANG_PROGRAM()], [break 2], [CXX=])
                    fi
                fi
            done
            AC_LANG_POP(C++)
            if test -n "$CXX"; then
                AC_MSG_RESULT($CXX)
            else
                AC_MSG_RESULT(failed)
                AC_MSG_FAILURE([can not determine the C++ compiler underlying $cit_MPICXX])
            fi
        fi
        cit_MPICXX_underlying_CXX=$CXX
    fi
fi
])dnl CIT_PROG_MPICXX
dnl end of file


# CIT_CHECK_LIB_MPI
# -----------------
AC_DEFUN([CIT_CHECK_LIB_MPI], [
# $Id: cit_mpi.m4 7980 2007-09-18 01:01:59Z leif $
AC_REQUIRE([_CIT_PROG_MPICC])dnl
AC_ARG_VAR(MPILIBS, [MPI linker flags, e.g. -L<mpi lib dir> -lmpi])
AC_SUBST(MPILIBS)
cit_mpi_save_CC=$CC
cit_mpi_save_CXX=$CXX
cit_mpi_save_LIBS=$LIBS
CC=$MPICC
CXX=$MPICXX
LIBS="$MPILIBS $LIBS"
# If MPILIBS is set, check to see if it works.
# If MPILIBS is not set, check to see if it is needed.
AC_CHECK_FUNC(MPI_Init, [], [
    if test -n "$MPILIBS"; then
        AC_MSG_ERROR([function MPI_Init not found; check MPILIBS])
    fi
    # MPILIBS is needed but was not set.
    AC_LANG_CASE(
        [C], [
            cit_mpicmd=$cit_MPICC
        ],
        [C++], [
            cit_mpicmd=$cit_MPICXX
            test -z "$cit_mpicmd" && cit_mpicmd=$cit_MPICC
        ]
    )
    if test -n "$cit_mpicmd"; then
        # Try to guess the correct value for MPILIBS using an MPI wrapper.
        CIT_MPI_LIBS(cit_libs, $cit_mpicmd, [
            LIBS="$cit_libs $cit_mpi_save_LIBS"
            unset ac_cv_func_MPI_Init
            AC_CHECK_FUNC(MPI_Init, [
                MPILIBS=$cit_libs
                export MPILIBS
            ], [
                _CIT_CHECK_LIB_MPI_FAILED
            ])
        ], [
            _CIT_CHECK_LIB_MPI_FAILED
        ])
    else
        # Desperate, last-ditch effort.
        cit_libs=
        for cit_lib in mpi mpich; do
            AC_CHECK_LIB($cit_lib, MPI_Init, [
                cit_libs="-l$cit_lib"
                MPILIBS=$cit_libs
                export MPILIBS
                break])
        done
        if test -z "$cit_libs"; then
            _CIT_CHECK_LIB_MPI_FAILED
        fi
    fi
])
LIBS=$cit_mpi_save_LIBS
CXX=$cit_mpi_save_CXX
CC=$cit_mpi_save_CC
])dnl CIT_CHECK_LIB_MPI


# _CIT_CHECK_LIB_MPI_FAILED
# -------------------------
AC_DEFUN([_CIT_CHECK_LIB_MPI_FAILED], [
AC_MSG_ERROR([no MPI library found

    Set the MPICC, MPICXX, MPIINCLUDES, and MPILIBS environment variables
    to specify how to build MPI programs.
])
])dnl _CIT_CHECK_LIB_MPI_FAILED


# CIT_HEADER_MPI
# --------------
AC_DEFUN([CIT_HEADER_MPI], [
# $Id: cit_mpi.m4 7980 2007-09-18 01:01:59Z leif $
AC_LANG_CASE(
    [C], [
        AC_CHECK_HEADER([mpi.h], [], [AC_MSG_ERROR([header 'mpi.h' not found])])
    ],
    [C++], [
        CIT_MPI_CHECK_CXX_LINK(cit_MPI_CPPFLAGS, [],
                               _CIT_TRIVIAL_MPI_PROGRAM,
                               [whether we can link a trivial C++ MPI program],
                               [],
                               AC_MSG_FAILURE([cannot link a trivial C++ MPI program using $CXX]))
        CPPFLAGS="$cit_MPI_CPPFLAGS $CPPFLAGS"
])
])dnl CIT_HEADER_MPI


# CIT_MPI_CHECK_CXX_LINK(INCLUDES, LIBS, PROGRAM,
#                        MSG, IF-WORKS, IF-NOT)
# -----------------------------------------------
AC_DEFUN([CIT_MPI_CHECK_CXX_LINK], [
# $Id: cit_mpi.m4 7980 2007-09-18 01:01:59Z leif $
AC_LANG_ASSERT(C++)
AC_MSG_CHECKING($4)
CIT_MPI_CXX_LINK_IFELSE(cit_arg, $$1, $2, $3,
[
    if test -z "$cit_arg"; then
	AC_MSG_RESULT(yes)
    else
	AC_MSG_RESULT([yes, with $cit_arg])
    fi
    $1="$cit_arg [$]$1"
    $5
], [
    AC_MSG_RESULT(no)
    $6
])
])


# CIT_MPI_CXX_LINK_IFELSE(DEFINES,
#                         INCLUDES, LIBS, PROGRAM,
#                         IF-WORKS, IF-NOT)
# ------------------------------------------------
# Verify that the MPI library is link-compatible with CXX (which could
# be different than the C++ compiler used to build the MPI library) by
# attempting to compile and link PROGRAM.  If there is a problem,
# attempt to work-around it by preventing MPI's C++ bindings from
# being #included.  If successful, set DEFINES to the preprocessor
# flags (if any) needed to successfully compile and link PROGRAM and
# evaluate IF-WORKS; otherwise, evaluate IF-NOT.
AC_DEFUN([CIT_MPI_CXX_LINK_IFELSE], [
# $Id: cit_mpi.m4 7980 2007-09-18 01:01:59Z leif $
AC_LANG_ASSERT(C++)
$1=
cit_mpi_cxx_link_save_LIBS=$LIBS
cit_mpi_cxx_link_save_CPPFLAGS=$CPPFLAGS
LIBS="$3 $LIBS"
CPPFLAGS="$2 $cit_mpi_cxx_link_save_CPPFLAGS"
AC_LINK_IFELSE([$4], [$5], [
    for cit_skip_mpicxx_define in CIT_SKIP_MPICXX_DEFINES
    do
	CPPFLAGS="$cit_skip_mpicxx_define $2 $cit_mpi_cxx_link_save_CPPFLAGS"
	AC_LINK_IFELSE([$4], [
	    $1=$cit_skip_mpicxx_define
            $5
	    break
	], [
            $6
	])
    done
])
CPPFLAGS=$cit_mpi_cxx_link_save_CPPFLAGS
LIBS=$cit_mpi_cxx_link_save_LIBS
])dnl CIT_MPI_CXX_LINK_IFELSE


# _CIT_TRIVIAL_MPI_PROGRAM
# ------------------------
AC_DEFUN([_CIT_TRIVIAL_MPI_PROGRAM], [
AC_LANG_PROGRAM([[
#include <stdio.h>
#include <mpi.h>
]], [[
    MPI_Init(0,0);
    MPI_Finalize();
]])
])dnl _CIT_TRIVIAL_MPI_PROGRAM


# CIT_MPI_LIBS(LIBS, COMMAND,
#              ACTION-IF-FOUND, ACTION-IF-NOT-FOUND)
# -------------------------------------------------------
# Guess the libraries used by the MPI wrapper.
AC_DEFUN([CIT_MPI_LIBS], [
# $Id: cit_mpi.m4 7980 2007-09-18 01:01:59Z leif $
AC_MSG_CHECKING([for the libraries used by $2])
$1=
for cit_arg_show in CIT_MPI_LINK_INFO_SWITCHES
do
    cit_cmd="$2 $cit_arg_show"
    if $cit_cmd >/dev/null 2>&1; then
	cit_args=`$cit_cmd 2>/dev/null`
	test -z "$cit_args" && continue
	for cit_arg in $cit_args
	do
	    case $cit_arg in
		-L* | -l* | -pthread* [)] $1="[$]$1 $cit_arg" ;;
	    esac
	done
	test -z "[$]$1" && continue
	break
    fi
done
if test -n "[$]$1"; then
    AC_MSG_RESULT([[$]$1])
    $3
else
    AC_MSG_RESULT(failed)
    $4
fi
])dnl CIT_MPI_LIBS


# CIT_MPI_INCLUDES(INCLUDES, COMMAND,
#                  ACTION-IF-FOUND, ACTION-IF-NOT-FOUND)
# ----------------------------------------------------------
# Guess the includes used by the MPI wrapper.
AC_DEFUN([CIT_MPI_INCLUDES], [
# $Id: cit_mpi.m4 7980 2007-09-18 01:01:59Z leif $
AC_MSG_CHECKING([for the includes used by $2])
$1=
for cit_arg_show in CIT_MPI_COMPILE_INFO_SWITCHES
do
    cit_cmd="$2 -c $cit_arg_show"
    if $cit_cmd >/dev/null 2>&1; then
	cit_args=`$cit_cmd 2>/dev/null`
	test -z "$cit_args" && continue
	for cit_arg in $cit_args
	do
	    case $cit_arg in
		-I* [)] $1="[$]$1 $cit_arg" ;;
	    esac
	done
	test -z "[$]$1" && continue
	break
    fi
done
if test -n "[$]$1"; then
    AC_MSG_RESULT([[$]$1])
    $3
else
    AC_MSG_RESULT(failed)
    $4
fi
])dnl CIT_MPI_INCLUDES


# CIT_MPI_COMPILE_INFO_SWITCHES
# CIT_MPI_LINK_INFO_SWITCHES
# -----------------------------
# The variety of flags used by MPICH, LAM/MPI, Open MPI, and ChaMPIon/Pro.
# NYI: mpxlc/mpcc (xlc?), mpcc_r (xlc_r?)
AC_DEFUN([CIT_MPI_COMPILE_INFO_SWITCHES], ["-show" "-showme" "-echo" "-compile_info"])
AC_DEFUN([CIT_MPI_LINK_INFO_SWITCHES], ["-show" "-showme" "-echo" "-link_info"])


# CIT_SKIP_MPICXX_DEFINES
# -----------------------
# Switches to disable inclusion of C++ MPI bindings.
AC_DEFUN([CIT_SKIP_MPICXX_DEFINES], ["-DMPICH_SKIP_MPICXX" "-UHAVE_MPI_CPP" "-DLAM_WANT_MPI2CPP=0" "-DLAM_BUILDING=1" "-DOMPI_WANT_CXX_BINDINGS=0" "-DOMPI_BUILDING=1"])


dnl end of file
