# -*- Autoconf -*-


## ---------------------------- ##
## Autoconf macros for Fortran. ##
## ---------------------------- ##


# _CIT_FC_MAIN
# ------------
# Define {F77,FC}_MAIN to the name of the alternate main() function
# for use with the Fortran libraries (i.e., MAIN__ or whatever), or
# 'main' if no such alternate name is found.
#
# As of Autoconf 2.59, the macro AC_FC_MAIN does not work with ifort
# v9, because the macro assumes that 'main' will be resolved by
# FCLIBS, but FCLIBS does not include Intel's 'for_main.o'.  This
# macro simply links with the Fortran compiler instead.
#
AC_DEFUN([_CIT_FC_MAIN],
[_AC_FORTRAN_ASSERT()dnl
AC_CACHE_CHECK([for alternate main to link with Fortran libraries],
               ac_cv_[]_AC_LANG_ABBREV[]_main,
[ac_[]_AC_LANG_ABBREV[]_m_save_LIBS=$LIBS
 LIBS="cfortran_test.$ac_objext $LIBS"
 ac_fortran_dm_var=[]_AC_FC[]_DUMMY_MAIN
 ac_cv_fortran_main="main" # default entry point name
 for ac_func in MAIN__ MAIN_ __main MAIN _MAIN __MAIN main_ main__ _main; do
   AC_LANG_PUSH(C)
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([@%:@ifdef FC_DUMMY_MAIN_EQ_F77
@%:@  undef F77_DUMMY_MAIN
@%:@  undef FC_DUMMY_MAIN
@%:@else
@%:@  undef $ac_fortran_dm_var
@%:@endif
@%:@define main $ac_func])],
                  [mv conftest.$ac_objext cfortran_test.$ac_objext],
                  [AC_MSG_FAILURE([cannot compile a simple C program])])
   AC_LANG_POP(C)
   AC_LINK_IFELSE([AC_LANG_SOURCE(
[      subroutine foobar()
      return
      end])], [ac_cv_fortran_main=$ac_func; break])
   rm -f cfortran_test* conftest*
 done
 ac_cv_[]_AC_LANG_ABBREV[]_main=$ac_cv_fortran_main
 rm -f cfortran_test* conftest*
 LIBS=$ac_[]_AC_LANG_ABBREV[]_m_save_LIBS
])
AC_DEFINE_UNQUOTED([]_AC_FC[]_MAIN, $ac_cv_[]_AC_LANG_ABBREV[]_main,
                   [Define to alternate name for `main' routine that is
                    called from a `main' in the Fortran libraries.])
])# _CIT_FC_MAIN


# CIT_F77_MAIN
# ------------
AC_DEFUN([CIT_F77_MAIN],
[AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])dnl
AC_LANG_PUSH(Fortran 77)dnl
_AC_FC_MAIN
AC_LANG_POP(Fortran 77)dnl
])# CIT_F77_MAIN


# CIT_FC_MAIN
# -----------
AC_DEFUN([CIT_FC_MAIN],
[AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])dnl
AC_LANG_PUSH(Fortran)dnl
_CIT_FC_MAIN
AC_LANG_POP(Fortran)dnl
])# CIT_FC_MAIN


# CIT_FC_OPEN_APPEND
# ------------------
AC_DEFUN([CIT_FC_OPEN_APPEND], [
AC_LANG_PUSH(Fortran)
cit_fc_append=no
AC_MSG_CHECKING([whether $FC supports OPEN control item 'position="append"'])
AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([], [[      open(10,file="foo",status="old",position="append")]])
], [
    AC_MSG_RESULT(yes)
    FCFLAGS="-DFORTRAN_POSITION_APPEND $FCFLAGS"; export FCFLAGS
    cit_fc_append=yes
], [
    AC_MSG_RESULT(no)
])
AC_MSG_CHECKING([whether $FC supports OPEN control item 'access="append"'])
AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([], [[      open(10,file="foo",status="old",access="append")]])
], [
    AC_MSG_RESULT(yes)
    FCFLAGS="-DFORTRAN_ACCESS_APPEND $FCFLAGS"; export FCFLAGS
    cit_fc_append=yes
], [
    AC_MSG_RESULT(no)
])
AS_IF([test $cit_fc_append = yes], [], [
    AC_MSG_FAILURE([cannot determine method for appending to Fortran files])
])
AC_LANG_POP(Fortran)
])dnl CIT_FC_OPEN_APPEND


# CIT_FC_STREAM_IO
# ----------------
AC_DEFUN([CIT_FC_STREAM_IO], [
AC_LANG_PUSH(Fortran)
AC_MSG_CHECKING([whether $FC supports stream i/o])
AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([], [[      open(10,file="foo",status="new",access="stream",
     & form="unformatted")
      write(10,pos=1) 1,2,3.0d0]])
], [
    AC_MSG_RESULT(yes)
    FCFLAGS="-DFORTRAN_STREAM_IO $FCFLAGS"; export FCFLAGS
], [
        AC_MSG_RESULT(no)
        AC_MSG_CHECKING([whether $FC supports f77-style binary direct-access i/o])
        AC_COMPILE_IFELSE([
            AC_LANG_PROGRAM([], [[      open(10,file="foo",status="new",access="direct",recl=1,
     & form="unformatted")
      write(10,rec=1) 1,2,3.0d0]])
    ], [
        AC_MSG_RESULT(yes)
        FCFLAGS="-DFORTRAN_F77_IO $FCFLAGS"; export FCFLAGS
        AC_MSG_CHECKING([whether $FC supports I/O specifiers 'advance' and 'eor'])
        AC_COMPILE_IFELSE([
            AC_LANG_PROGRAM([], [[      open(10,file="foo",status="new",access="direct",recl=1,
     & form="unformatted")
      write(10,rec=1,advance='yes',eor=10) 1,2,3.0d0
 10   continue]])
        ], [
            AC_MSG_RESULT(yes)
            FCFLAGS="-DFORTRAN_EOR $FCFLAGS"; export FCFLAGS
        ], [
            AC_MSG_RESULT(no)
        ])
    ], [
        AC_MSG_RESULT(no)
        AC_MSG_WARN([cannot determine how to produce binary direct-access files with variable record length])
        FCFLAGS="-DFORTRAN_NO_BINARY $FCFLAGS"; export FCFLAGS
    ])
])
AC_LANG_POP(Fortran)
])dnl CIT_FC_STREAM_IO


# CIT_FC_MPI_MODULE(FILENAME, MPIFC, MPIFCFLAGS)
# -----------------------------------------------------
AC_DEFUN([CIT_FC_MPI_MODULE], [
# Use 'mpi' module or 'mpif.h', as appropriate.  UNFINISHED.  This
# strategy doesn't play well with "implicit none": whether the
# generated header must be included before or after "implicit none"
# depends upon the result of the test!  It might be possible to make
# "use mpi" always work: simply generate an 'mpi' module if the MPI
# library doesn't provide one.  The generated module would simply
# "include 'mpif.h'".
AC_LANG_PUSH(Fortran)

ofile=$1
cfgfile="${ofile}T"
trap "rm \"$cfgfile\"; exit 1" 1 2 15
rm -f "$cfgfile"

cit_fc_save_fc=$FC
cit_fc_save_fcflags=$FCFLAGS
FC=$2
FCFLAGS="$FCFLAGS $3"

AC_MSG_CHECKING([whether "use mpi" works])
AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([], [[
      use mpi
      integer ier
      call MPI_INIT(ier)
      call MPI_FINALIZE(ier)
]])
], [
    AC_MSG_RESULT(yes)
    cit_fc_header="use mpi"
], [
    AC_MSG_RESULT(no)
    AC_MSG_CHECKING([whether mpif.h works])
    AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM([], [[
      include 'mpif.h'
      integer ier
      call MPI_INIT(ier)
      call MPI_FINALIZE(ier)
]])
    ], [
        AC_MSG_RESULT(yes)
dnl Allow projects to simply include the standard 'mpif.h' everywhere.
dnl If FILENAME is 'mpif.h', this macro will conditionally create a header
dnl to override the system header.
        if test "$ofile" = "mpif.h"; then
            cit_fc_header=none
        else
            cit_fc_header="include 'mpif.h'"
        fi
    ], [
        AC_MSG_RESULT(no)
        AC_MSG_FAILURE([cannot compile a trivial MPI program using $2])
    ])
])

if test "$cit_fc_header" != "none"; then
    AC_MSG_NOTICE([creating $ofile])
    cat >"$cfgfile" <<END_OF_HEADER
! $ofile.  Generated by configure.

      $cit_fc_header

END_OF_HEADER
    mv -f "$cfgfile" "$ofile" || \
        (rm -f "$ofile" && cp "$cfgfile" "$ofile" && rm -f "$cfgfile")
fi


FC=$cit_fc_save_fc
FCFLAGS=$cit_fc_save_fcflags

AC_LANG_POP(Fortran)
])dnl CIT_FC_MPI_MODULE


# CIT_FC_MPI_HEADER(MPIFC, MPIFCFLAGS)
# -----------------------------------------------------
AC_DEFUN([CIT_FC_MPI_HEADER], [
# Generate a Fortran 9x-compatible 'mpif.h', if necessary.
AC_LANG_PUSH(Fortran)

ofile="mpif.h"
cfgfile="${ofile}T"
trap "rm \"$cfgfile\"; exit 1" 1 2 15
rm -f "$cfgfile"

cit_fc_save_fc=$FC
cit_fc_save_fcflags=$FCFLAGS
FC=$1
FCFLAGS="$FCFLAGS $2"

AC_MSG_CHECKING([whether mpif.h works])
AC_COMPILE_IFELSE(_CIT_FC_TRIVIAL_MPI_PROGRAM, [
    AC_MSG_RESULT(yes)
], [
    AC_MSG_RESULT(no)
    cit_mpif_h=unknown
    cit_mpifc_info=`$FC -compile_info 2>/dev/null`
    for cit_arg in $cit_mpifc_info; do
        case $cit_arg in
            */mpif.h) cit_mpif_h="$cit_arg"; break;;
        esac
    done
    if test "$cit_mpif_h" == "unknown"; then
        AC_MSG_FAILURE([cannot compile a trivial MPI program using $1])
    fi

dnl Special hack for MPICH.
    AC_MSG_NOTICE([creating $ofile])
    cat >"$cfgfile" <<END_OF_HEADER
! $ofile.  Generated from $cit_mpif_h by configure.

END_OF_HEADER
    grep -v MPI_DISPLACEMENT_CURRENT "$cit_mpif_h" >>"$cfgfile"
    mv -f "$cfgfile" "$ofile" || \
        (rm -f "$ofile" && cp "$cfgfile" "$ofile" && rm -f "$cfgfile")

    AC_MSG_CHECKING([whether generated mpif.h works])
    AC_COMPILE_IFELSE(_CIT_FC_TRIVIAL_MPI_PROGRAM, [
        AC_MSG_RESULT(yes)
    ], [
        AC_MSG_RESULT(no)
        AC_MSG_FAILURE([cannot compile a trivial MPI program using $1])
    ])

])

FC=$cit_fc_save_fc
FCFLAGS=$cit_fc_save_fcflags

AC_LANG_POP(Fortran)
])dnl CIT_FC_MPI_HEADER


# _CIT_FC_TRIVIAL_MPI_PROGRAM
# ------------------------
AC_DEFUN([_CIT_FC_TRIVIAL_MPI_PROGRAM], [
AC_LANG_PROGRAM([], [[
      include 'mpif.h'
      integer, parameter :: CUSTOM_MPI_TYPE = MPI_REAL
      integer ier
      call MPI_INIT(ier)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
      call MPI_FINALIZE(ier)
]])
])dnl _CIT_FC_TRIVIAL_MPI_PROGRAM


dnl end of file
