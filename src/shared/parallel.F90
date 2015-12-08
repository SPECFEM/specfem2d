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

! parallel routines

!-------------------------------------------------------------------------------------------------
!
! MPI wrapper functions
!
!-------------------------------------------------------------------------------------------------

  subroutine init_mpi()

#ifdef USE_MPI
! standard include of the MPI library
  use mpi
#endif

  implicit none

#ifdef USE_MPI
  ! local parameters
  integer :: sizeprocs,myrank
  integer :: ier

  ! parallel version
  call MPI_INIT(ier)
  if (ier /= 0 ) stop 'Error initializing MPI'

  ! checks if getting size works
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  if (ier /= 0 ) stop 'Error getting MPI size'

  ! checks if getting rank works
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)
  if (ier /= 0 ) stop 'Error getting MPI rank'
#endif

  end subroutine init_mpi


!
!-------------------------------------------------------------------------------------------------
!

  subroutine finalize_mpi()

#ifdef USE_MPI
! standard include of the MPI library
  use mpi
#endif

  implicit none

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  ! synchronizes all
  call MPI_BARRIER(MPI_COMM_WORLD,ier)

  ! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)
  if (ier /= 0) stop 'Error finalizing MPI'
#endif

  end subroutine finalize_mpi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine abort_mpi()

#ifdef USE_MPI
  use mpi
#endif

  implicit none

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  ! stop all the MPI processes, and exit
  ! note: MPI_ABORT does not return, it makes the program exit with an error code of 30
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)
#endif

  stop 'error, program ended in exit_MPI'

  end subroutine abort_mpi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine synchronize_all()

#ifdef USE_MPI
! standard include of the MPI library
  use mpi
#endif

  implicit none

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_BARRIER(MPI_COMM_WORLD,ier)
  if (ier /= 0 ) stop 'Error synchronize MPI processes'
#endif

  end subroutine synchronize_all

!
!-------------------------------------------------------------------------------------------------
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


!-------------------------------------------------------------------------------------------------
!
! Send/Receive MPI
!
!-------------------------------------------------------------------------------------------------


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
!-------------------------------------------------------------------------------------------------
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


!-------------------------------------------------------------------------------------------------
!
! MPI math helper
!
!-------------------------------------------------------------------------------------------------


  subroutine min_all_i(sendbuf, recvbuf)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer:: sendbuf, recvbuf

#ifdef USE_MPI
  ! local parameters
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_MIN,0,MPI_COMM_WORLD,ier)
#else
  recvbuf = sendbuf
#endif

  end subroutine min_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine min_all_all_dp(sendbuf, recvbuf)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  double precision :: sendbuf, recvbuf

#ifdef USE_MPI
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ier)
#else
  recvbuf = sendbuf
#endif

  end subroutine min_all_all_dp

!
!-------------------------------------------------------------------------------------------------
!


  subroutine max_all_i(sendbuf, recvbuf)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer :: sendbuf, recvbuf

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ier)
#else
  recvbuf = sendbuf
#endif

  end subroutine max_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_all_dp(sendbuf, recvbuf)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  double precision :: sendbuf, recvbuf

#ifdef USE_MPI
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ier)
#else
  recvbuf = sendbuf
#endif

  end subroutine max_all_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine any_all_l(sendbuf, recvbuf)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  logical :: sendbuf, recvbuf

#ifdef USE_MPI
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ier)
#else
  recvbuf = sendbuf
#endif

  end subroutine any_all_l

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_dp(sendbuf, recvbuf)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  double precision :: sendbuf, recvbuf

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ier)
#else
    recvbuf = sendbuf
#endif

  end subroutine sum_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_all_dp(sendbuf, recvbuf)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  double precision :: sendbuf, recvbuf

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ier)
#else
    recvbuf = sendbuf
#endif

  end subroutine sum_all_all_dp


!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_cr(sendbuf, recvbuf)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  include "constants.h"
#ifdef USE_MPI
  include "precision.h"
#endif

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
#else
    recvbuf = sendbuf
#endif

  end subroutine sum_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_i(sendbuf, recvbuf)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer :: sendbuf, recvbuf

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ier)
#else
    recvbuf = sendbuf
#endif

  end subroutine sum_all_i

!
!-------------------------------------------------------------------------------------------------
!


  subroutine sum_all_all_i(sendbuf, recvbuf)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer :: sendbuf, recvbuf

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ier)
#else
  recvbuf = sendbuf
#endif

  end subroutine sum_all_all_i

!-------------------------------------------------------------------------------------------------
!
! MPI world helper
!
!-------------------------------------------------------------------------------------------------

  subroutine world_size(sizeval)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer,intent(out) :: sizeval

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeval,ier)
  if (ier /= 0 ) stop 'Error getting MPI world size'
#else
  ! single process
  sizeval = 1
#endif

  end subroutine world_size

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_rank(rank)

#ifdef USE_MPI
! standard include of the MPI library
  use mpi
#endif

  implicit none

  integer,intent(out) :: rank

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ier)
  if (ier /= 0 ) stop 'Error getting MPI rank'
#else
  ! always returns master rank zero
  rank = 0
#endif

  end subroutine world_rank
