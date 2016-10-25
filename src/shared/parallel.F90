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

  stop 'Error, program ended in exit_MPI'

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


  subroutine wait_req(req)

#ifdef USE_MPI
! standard include of the MPI library
  use mpi
#endif

  implicit none

  integer :: req

#ifndef USE_MPI
  integer :: dummy
#endif

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_WAIT(req,MPI_STATUS_IGNORE,ier)
#else
  ! to avoid compiler warning
  dummy = req
#endif

  end subroutine wait_req


!-------------------------------------------------------------------------------------------------
!
! MPI broadcasting helper
!
!-------------------------------------------------------------------------------------------------


  subroutine bcast_all_i(buffer, countval)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer :: countval
  integer, dimension(countval) :: buffer
#ifndef USE_MPI
  integer :: dummy
#endif

#ifdef USE_MPI
  integer :: ier

  call MPI_BCAST(buffer,countval,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
#else
  ! to avoid compiler warning
  dummy = buffer(1)
#endif

  end subroutine bcast_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singlei(buffer)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer :: buffer
#ifndef USE_MPI
  integer :: dummy
#endif

#ifdef USE_MPI
  integer :: ier

  call MPI_BCAST(buffer,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
#else
  ! to avoid compiler warning
  dummy = buffer
#endif

  end subroutine bcast_all_singlei

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singlel(buffer)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  logical :: buffer
#ifndef USE_MPI
  logical :: dummy
#endif

#ifdef USE_MPI
  integer :: ier

  call MPI_BCAST(buffer,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
#else
  ! to avoid compiler warning
  dummy = buffer
#endif

  end subroutine bcast_all_singlel


!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singledp(buffer)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  double precision :: buffer
#ifndef USE_MPI
  double precision :: dummy
#endif

#ifdef USE_MPI
  integer :: ier

  call MPI_BCAST(buffer,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
#else
  ! to avoid compiler warning
  dummy = buffer
#endif

  end subroutine bcast_all_singledp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_string(buffer)

#ifdef USE_MPI
  use mpi
#endif

  use constants, only: MAX_STRING_LEN

  implicit none

  character(len=MAX_STRING_LEN) :: buffer

#ifndef USE_MPI
  character(len=MAX_STRING_LEN) :: dummy
#endif

#ifdef USE_MPI
  integer :: ier

  call MPI_BCAST(buffer,MAX_STRING_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
#else
  ! to avoid compiler warning
  dummy = buffer
#endif

  end subroutine bcast_all_string


!-------------------------------------------------------------------------------------------------
!
! Send/Receive MPI
!
!-------------------------------------------------------------------------------------------------


  subroutine isend_cr(sendbuf, sendcount, dest, sendtag, req)

! asynchronuous send

#ifdef USE_MPI
! standard include of the MPI library
  use mpi
#endif

  use constants, only: CUSTOM_REAL

  implicit none

#ifdef USE_MPI
  include "precision.h"
#endif

  integer :: sendcount, dest, sendtag, req
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf

#ifndef USE_MPI
  integer :: dummy
#endif

#ifdef USE_MPI
  integer :: ier

  call MPI_ISEND(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag,MPI_COMM_WORLD,req,ier)
#else
  stop 'isend_cr not implemented for serial code'
  ! to avoid compiler warning
  dummy = sendbuf(1)
  dummy = dest
  dummy = sendtag
  dummy = req
#endif

  end subroutine isend_cr


!
!-------------------------------------------------------------------------------------------------
!
#ifdef USE_MPI

  subroutine send_singlei(sendbuf, dest, sendtag)

! synchronuous/blocking send

  use mpi

  implicit none

  integer :: dest,sendtag
  integer :: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf,1,MPI_INTEGER,dest,sendtag,MPI_COMM_WORLD,ier)

  end subroutine send_singlei

#endif

!
!-------------------------------------------------------------------------------------------------
!
#ifdef USE_MPI

  subroutine send_singledp(sendbuf, dest, sendtag)

! synchronuous/blocking send

  use mpi

  implicit none

  integer :: dest,sendtag
  double precision :: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf,1,MPI_DOUBLE_PRECISION,dest,sendtag,MPI_COMM_WORLD,ier)

  end subroutine send_singledp

#endif

!
!-------------------------------------------------------------------------------------------------
!

  subroutine irecv_cr(recvbuf, recvcount, dest, recvtag, req)

! asynchronuous receive

#ifdef USE_MPI
! standard include of the MPI library
  use mpi
#endif

  use constants, only: CUSTOM_REAL

  implicit none

#ifdef USE_MPI
  include "precision.h"
#endif

  integer :: recvcount, dest, recvtag, req
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf
#ifndef USE_MPI
  integer :: dummy
#endif

#ifdef USE_MPI
  integer :: ier

  call MPI_IRECV(recvbuf,recvcount,CUSTOM_MPI_TYPE,dest,recvtag,MPI_COMM_WORLD,req,ier)
#else
  stop 'irecv_cr not implemented for serial code'
  ! to avoid compiler warning
  dummy = recvbuf(1)
  dummy = dest
  dummy = recvtag
  dummy = req
#endif

  end subroutine irecv_cr

!
!-------------------------------------------------------------------------------------------------
!
#ifdef USE_MPI

  subroutine recv_singlei(recvbuf, dest, recvtag )

! synchronuous/blocking receive

  use mpi

  implicit none

  integer :: dest,recvtag
  integer :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,1,MPI_INTEGER,dest,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

  end subroutine recv_singlei

#endif

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_MPI

  subroutine recv_singledp(recvbuf, dest, recvtag )

! synchronuous/blocking receive

  use mpi

  implicit none

  integer :: dest,recvtag
  double precision :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,1,MPI_DOUBLE_PRECISION,dest,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

  end subroutine recv_singledp

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

  subroutine max_all_cr(sendbuf, recvbuf)

#ifdef USE_MPI
  use mpi
#endif
  use constants, only: CUSTOM_REAL

  implicit none

#ifdef USE_MPI
  include "precision.h"
#endif

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf

#ifdef USE_MPI
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
#else
  recvbuf = sendbuf
#endif

  end subroutine max_all_cr

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
  use constants, only: CUSTOM_REAL

  implicit none

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
! MPI gather helper
!
!-------------------------------------------------------------------------------------------------

  subroutine gather_all_i(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer :: sendcnt, recvcount, NPROC
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcount,0:NPROC-1) :: recvbuf

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_GATHER(sendbuf,sendcnt,MPI_INTEGER,recvbuf,recvcount,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
#else
  recvbuf(:,0) = sendbuf(:)
#endif

  end subroutine gather_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_dp(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer :: sendcnt, recvcount, NPROC
  double precision, dimension(sendcnt) :: sendbuf
  double precision, dimension(recvcount,0:NPROC-1) :: recvbuf

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_GATHER(sendbuf,sendcnt,MPI_DOUBLE_PRECISION,recvbuf,recvcount,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
#else
  recvbuf(:,0) = sendbuf(:)
#endif

  end subroutine gather_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_singlei(sendbuf, recvbuf, NPROC)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer :: NPROC
  integer :: sendbuf
  integer, dimension(0:NPROC-1) :: recvbuf

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_GATHER(sendbuf,1,MPI_INTEGER,recvbuf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
#else
  recvbuf(0) = sendbuf
#endif

  end subroutine gather_all_singlei

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_all_singlei(sendbuf, recvbuf, NPROC)

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer :: NPROC
  integer :: sendbuf
  integer, dimension(0:NPROC-1) :: recvbuf

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_ALLGATHER(sendbuf,1,MPI_INTEGER,recvbuf,1,MPI_INTEGER,MPI_COMM_WORLD,ier)
#else
  recvbuf(0) = sendbuf
#endif

  end subroutine gather_all_all_singlei


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
