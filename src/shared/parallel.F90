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
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and Inria at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

! parallel routines

!-------------------------------------------------------------------------------------------------
!
! MPI wrapper functions
!
!-------------------------------------------------------------------------------------------------

  subroutine init_mpi(NPROC,myrank)

#ifdef USE_MPI
! standard include of the MPI library
  use mpi
#endif

  implicit none

  integer,intent(out) :: NPROC,myrank

#ifdef USE_MPI
  integer :: ier

  ! parallel version
  call MPI_INIT(ier)
  if (ier /= 0 ) call exit_MPI('Error MPI initialization')

  call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,ier)
  if (ier /= 0 ) call exit_MPI('Error getting MPI size')

  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)
  if (ier /= 0 ) call exit_MPI('Error getting MPI rank')

#else

  ! serial version
  ! compilation without MPI support -DUSE_MPI
  NPROC = 1
  myrank = 0

#endif

  end subroutine init_mpi


!
!----
!

  subroutine finalize_mpi()

#ifdef USE_MPI
! standard include of the MPI library
  use mpi
#endif

  implicit none

#ifdef USE_MPI
  integer :: ier

  ! synchronizes all
  call MPI_BARRIER(MPI_COMM_WORLD,ier)

  ! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)
  if (ier /= 0) stop 'Error finalizing MPI'
#endif

  end subroutine finalize_mpi

!
!----
!

  subroutine synchronize_all()

#ifdef USE_MPI
! standard include of the MPI library
  use mpi
#endif

  implicit none

#ifdef USE_MPI
  integer ier

  call MPI_BARRIER(MPI_COMM_WORLD,ier)
  if (ier /= 0 ) stop 'Error synchronize MPI processes'
#endif

  end subroutine synchronize_all

!
!-------------------------------------------------------------------------------------------------
!


#ifdef USE_MPI

  subroutine isend_cr(sendbuf, sendcount, dest, sendtag, req)

! standard include of the MPI library
  use mpi

  implicit none

  include "constants.h"
  include "precision.h"

  integer sendcount, dest, sendtag, req
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf

  integer ier

  call MPI_ISEND(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag,MPI_COMM_WORLD,req,ier)

  end subroutine isend_cr

#endif

!
!----
!

#ifdef USE_MPI

  subroutine irecv_cr(recvbuf, recvcount, dest, recvtag, req)

! standard include of the MPI library
  use mpi

  implicit none

  include "constants.h"
  include "precision.h"

  integer recvcount, dest, recvtag, req
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

  integer ier

  call MPI_IRECV(recvbuf,recvcount,CUSTOM_MPI_TYPE,dest,recvtag,MPI_COMM_WORLD,req,ier)

  end subroutine irecv_cr

#endif

!
!----
!

#ifdef USE_MPI

  subroutine wait_req(req)

! standard include of the MPI library
  use mpi

  implicit none

  integer :: req

  integer, dimension(MPI_STATUS_SIZE) :: req_mpi_status

  integer :: ier

  call mpi_wait(req,req_mpi_status,ier)

  end subroutine wait_req

#endif


!
!----
!

#ifdef USE_MPI

  subroutine min_all_i(sendbuf, recvbuf)

! standard include of the MPI library
  use mpi

  implicit none

  include "constants.h"

  integer:: sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_MIN,0,MPI_COMM_WORLD,ier)

  end subroutine min_all_i

#endif

!
!----
!

#ifdef USE_MPI

  subroutine max_all_i(sendbuf, recvbuf)

! standard include of the MPI library
  use mpi

  implicit none

  include "constants.h"

  integer :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ier)

  end subroutine max_all_i

#endif

