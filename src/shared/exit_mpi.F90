!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

!-----------------------------------------------
! subroutine to stop the code whether sequential or parallel.
!-----------------------------------------------
  subroutine exit_MPI(myrank,error_msg)

#ifdef USE_MPI
  use mpi
#endif
  use constants, only: MAX_STRING_LEN,IMAIN,ISTANDARD_OUTPUT

  implicit none

  ! identifier for error message file
  integer, parameter :: IERROR = 30

  integer :: myrank
  character(len=*) :: error_msg

  character(len=MAX_STRING_LEN) :: outputname

  ! write error message to screen
  write(*,*) error_msg(1:len(error_msg))
  write(*,*) 'Error detected, aborting MPI... proc ',myrank

! write error message to file
  write(outputname,"('/error_message',i6.6,'.txt')") myrank
  open(unit=IERROR,file='OUTPUT_FILES'//outputname,status='unknown')
  write(IERROR,*) error_msg(1:len(error_msg))
  write(IERROR,*) 'Error detected, aborting MPI... proc ',myrank
  close(IERROR)

! close output file
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) close(IMAIN)

  call abort_mpi()

  end subroutine exit_MPI

!-------------------------------------------------------------------------------------------------
!
! I/O wrapper function
!
!-------------------------------------------------------------------------------------------------

  subroutine flush_IMAIN()

  use constants, only: IMAIN

  implicit none

  ! only master process writes out to main output file
  ! file I/O in fortran is buffered by default
  !
  ! note: Fortran2003 includes a FLUSH statement
  !          which is implemented by most compilers by now
  !
  ! otherwise:
  !   a) comment out the line below
  !   b) try to use instead: call flush(IMAIN)

  flush(IMAIN)

  end subroutine flush_IMAIN

