!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

! We added the ability to run several calculations (several earthquakes)
! in an embarrassingly-parallel fashion from within the same run;
! this can be useful when using a very large supercomputer to compute
! many earthquakes in a catalog, in which case it can be better from
! a batch job submission point of view to start fewer and much larger jobs,
! each of them computing several earthquakes in parallel.
! To turn that option on, set parameter NUMBER_OF_SIMULTANEOUS_RUNS to a value greater than 1 in the Par_file.
! To implement that, we create NUMBER_OF_SIMULTANEOUS_RUNS MPI sub-communicators,
! each of them being labeled "my_local_mpi_comm_world", and we use them
! in all the routines in "src/shared/parallel.f90", except in MPI_ABORT() because in that case
! we need to kill the entire run.
! When that option is on, of course the number of processor cores used to start
! the code in the batch system must be a multiple of NUMBER_OF_SIMULTANEOUS_RUNS,
! all the individual runs must use the same number of processor cores,
! which as usual is NPROC in the input file DATA/Par_file,
! and thus the total number of processor cores to request from the batch system
! should be NUMBER_OF_SIMULTANEOUS_RUNS * NPROC.
! All the runs to perform must be placed in directories called run0001, run0002, run0003 and so on
! (with exactly four digits).

!-------------------------------------------------------------------------------------------------
!
! Parallel routines.  All MPI calls should belong in this file!
!
!-------------------------------------------------------------------------------------------------

module my_mpi_communicator

! main parameter module for specfem simulations

  !use MPI ! TODO

  implicit none

  integer :: my_local_mpi_comm_world, my_local_mpi_comm_for_bcast

end module my_mpi_communicator

!-------------------------------------------------------------------------------------------------
!
! MPI wrapper functions
!
!-------------------------------------------------------------------------------------------------

  subroutine init_mpi()

  use my_mpi_communicator
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

#ifdef USE_MPI
! standard include of the MPI library
  use mpi
#endif

  implicit none

#ifdef USE_MPI
  integer :: myrank
  ! local parameters
  integer :: sizeprocs
  integer :: ier

  ! parallel version
  call MPI_INIT(ier)
  if (ier /= 0 ) call stop_the_code('Error initializing MPI')

  ! checks if getting size works
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  if (ier /= 0 ) call stop_the_code('Error getting MPI size')
  ! we need to make sure that NUMBER_OF_SIMULTANEOUS_RUNS and BROADCAST_SAME_MESH_AND_MODEL are read before calling world_split()
  ! thus read the parameter file
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)
  if (ier /= 0 ) call stop_the_code('Error getting MPI rank')
  if (myrank == 0) then
    call open_parameter_file_from_master_only()
    ! we need to make sure that NUMBER_OF_SIMULTANEOUS_RUNS and BROADCAST_SAME_MESH_AND_MODEL are read
    call read_value_integer_p(NUMBER_OF_SIMULTANEOUS_RUNS, 'NUMBER_OF_SIMULTANEOUS_RUNS')
    call read_value_logical_p(BROADCAST_SAME_MESH_AND_MODEL, 'BROADCAST_SAME_MESH_AND_MODEL')
    ! close parameter file
    call close_parameter_file()
  endif
  ! broadcast parameters read from master to all processes
  my_local_mpi_comm_world = MPI_COMM_WORLD
  call bcast_all_singlei(NUMBER_OF_SIMULTANEOUS_RUNS)
  call bcast_all_singlel(BROADCAST_SAME_MESH_AND_MODEL)
! create sub-communicators if needed, if running more than one earthquake from the same job
  call world_split()
#else
  NUMBER_OF_SIMULTANEOUS_RUNS = NUMBER_OF_SIMULTANEOUS_RUNS ! To avoid compiler warning
  BROADCAST_SAME_MESH_AND_MODEL = BROADCAST_SAME_MESH_AND_MODEL ! To avoid compiler warning
  ! we need to make sure that NUMBER_OF_SIMULTANEOUS_RUNS is read, thus read the parameter file
  ! initialize
  call read_parameter_file_init()
  ! open the Par_file
  call open_parameter_file()
  ! read only parameters (without receiver-line section, material tables or region definitions)
  call read_parameter_file_only()
#endif

  end subroutine init_mpi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine finalize_mpi()

  use my_mpi_communicator
#ifdef USE_MPI
  ! standard include of the MPI library
  use mpi
#endif

  implicit none

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  ! close sub-communicators if needed, if running more than one earthquake from the same job
  call world_unsplit()

  ! synchronizes all
  call MPI_BARRIER(MPI_COMM_WORLD,ier)

  ! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)
  if (ier /= 0) call stop_the_code('Error finalizing MPI')
#endif

  end subroutine finalize_mpi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine abort_mpi()

  use my_mpi_communicator
#ifdef USE_MPI
  use mpi
  use constants, only: MAX_STRING_LEN,mygroup
  use shared_input_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS
#endif

  implicit none

#ifdef USE_MPI
  ! local parameters

  integer :: my_local_rank,my_global_rank,ier
  logical :: run_file_exists
  character(len=MAX_STRING_LEN) :: filename

  ! get my local rank and my global rank (in the case of simultaneous jobs, for which we split
  ! the MPI communicator, they will be different; otherwise they are the same)
  call world_rank(my_local_rank)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_global_rank,ier)

  ! write a stamp file to disk to let the user know that the run failed
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) then
    ! notifies which run directory failed
    write(filename,"('run',i4.4,'_failed')") mygroup + 1
    inquire(file=trim(filename), exist=run_file_exists)
    if (run_file_exists) then
      open(unit=9765,file=trim(filename),status='old',position='append',action='write',iostat=ier)
    else
      open(unit=9765,file=trim(filename),status='new',action='write',iostat=ier)
    endif
    if (ier == 0) then
      write(9765,*) 'run ',mygroup+1,' with local rank ',my_local_rank,' and global rank ',my_global_rank,' failed'
      close(9765)
    endif

    ! notifies which rank failed
    write(filename,"('run_with_local_rank_',i8.8,'and_global_rank_',i8.8,'_failed')") my_local_rank,my_global_rank
    open(unit=9765,file=trim(filename),status='unknown',action='write')
    write(9765,*) 'run with local rank ',my_local_rank,' and global rank ',my_global_rank,' failed'
    close(9765)
  endif

  ! note: MPI_ABORT does not return, it makes the program exit with an error code of 30
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)
#endif

  call stop_the_code('Error, program ended in exit_MPI')

  end subroutine abort_mpi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine synchronize_all()

  use my_mpi_communicator
#ifdef USE_MPI
! standard include of the MPI library
  use mpi
#endif

  implicit none

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_BARRIER(my_local_mpi_comm_world,ier)
  if (ier /= 0 ) call stop_the_code('Error synchronize MPI processes')
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

  use my_mpi_communicator
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

  call MPI_BCAST(buffer,countval,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
#else
  ! to avoid compiler warning
  dummy = buffer(1)
#endif

  end subroutine bcast_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singlei(buffer)

 use my_mpi_communicator
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

  call MPI_BCAST(buffer,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
#else
  ! to avoid compiler warning
  dummy = buffer
#endif

  end subroutine bcast_all_singlei

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singlel(buffer)

  use my_mpi_communicator
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

  call MPI_BCAST(buffer,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
#else
  ! to avoid compiler warning
  dummy = buffer
#endif

  end subroutine bcast_all_singlel


!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singledp(buffer)

  use my_mpi_communicator
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

  call MPI_BCAST(buffer,1,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_world,ier)
#else
  ! to avoid compiler warning
  dummy = buffer
#endif

  end subroutine bcast_all_singledp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_string(buffer)

  use my_mpi_communicator
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

  call MPI_BCAST(buffer,MAX_STRING_LEN,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
#else
  ! to avoid compiler warning
  dummy = buffer
#endif

  end subroutine bcast_all_string

!
!---- broadcast to send the mesh and model to other simultaneous runs
!

  subroutine bcast_all_i_for_database(buffer, countval)

  use my_mpi_communicator
#ifdef USE_MPI
  use mpi
#endif
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none
  integer :: countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  integer :: buffer
#ifdef USE_MPI
  integer ier
#endif

  if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

#ifdef USE_MPI
  call MPI_BCAST(buffer,countval,MPI_INTEGER,0,my_local_mpi_comm_for_bcast,ier)
#else
  ! to avoid compiler warning
  buffer = buffer
  countval = countval
#endif

  end subroutine bcast_all_i_for_database

!
!-------------------------------------------------------------------------------------------------
!
! Send/Receive MPI
!
!-------------------------------------------------------------------------------------------------

  subroutine bcast_all_l_for_database(buffer, countval)

  use my_mpi_communicator
#ifdef USE_MPI
  use mpi
#endif
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  logical :: buffer

#ifdef USE_MPI
  integer ier
#endif

  if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

#ifdef USE_MPI
  call MPI_BCAST(buffer,countval,MPI_INTEGER,0,my_local_mpi_comm_for_bcast,ier)
#else
  ! to avoid compiler warning
  buffer = buffer
  countval = countval
#endif

  end subroutine bcast_all_l_for_database

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_cr_for_database(buffer, countval)

  use my_mpi_communicator
#ifdef USE_MPI
  use mpi
#endif
  use constants, only: CUSTOM_REAL
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none
#ifdef USE_MPI
  include "precision.h"
#endif

  integer countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  real(kind=CUSTOM_REAL) :: buffer

#ifdef USE_MPI
  integer ier
#endif

  if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

#ifdef USE_MPI
  call MPI_BCAST(buffer,countval,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_for_bcast,ier)
#else
  ! to avoid compiler warning
  buffer = buffer
  countval = countval
#endif

  end subroutine bcast_all_cr_for_database

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_dp_for_database(buffer, countval)

  use my_mpi_communicator
#ifdef USE_MPI
  use mpi
#endif
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  double precision :: buffer
#ifdef USE_MP
  integer ier
#endif

  if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

#ifdef USE_MP
  call MPI_BCAST(buffer,countval,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_for_bcast,ier)
#else
  ! to avoid compiler warning
  buffer = buffer
  countval = countval
#endif

  end subroutine bcast_all_dp_for_database

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_r_for_database(buffer, countval)

  use my_mpi_communicator
#ifdef USE_MPI
  use mpi
#endif
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  real :: buffer

#ifdef USE_MPI
  integer ier
#endif

  if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

#ifdef USE_MPI
  call MPI_BCAST(buffer,countval,MPI_REAL,0,my_local_mpi_comm_for_bcast,ier)
#else
  ! to avoid compiler warning
  buffer = buffer
  countval = countval
#endif

  end subroutine bcast_all_r_for_database

!-------------------------------------------------------------------------------------------------
!
! MPI math helper
!
!-------------------------------------------------------------------------------------------------

  subroutine min_all_i(sendbuf, recvbuf)

  use my_mpi_communicator
#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer:: sendbuf, recvbuf

#ifdef USE_MPI
  ! local parameters
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_MIN,0,my_local_mpi_comm_world,ier)
#else
  recvbuf = sendbuf
#endif

  end subroutine min_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_i(sendbuf, recvbuf)

  use my_mpi_communicator
#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer :: sendbuf, recvbuf

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_MAX,0,my_local_mpi_comm_world,ier)
#else
  recvbuf = sendbuf
#endif

  end subroutine max_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine min_all_all_dp(sendbuf, recvbuf)

  use my_mpi_communicator
#ifdef USE_MPI
  use mpi
#endif

  implicit none

  double precision :: sendbuf, recvbuf

#ifdef USE_MPI
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MIN,my_local_mpi_comm_world,ier)
#else
  recvbuf = sendbuf
#endif

  end subroutine min_all_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_cr(sendbuf, recvbuf)

  use my_mpi_communicator
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

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_MAX,0,my_local_mpi_comm_world,ier)
#else
  recvbuf = sendbuf
#endif

  end subroutine max_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_all_i(sendbuf, recvbuf)

  use my_mpi_communicator
#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer :: sendbuf, recvbuf

#ifdef USE_MPI
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_MAX,my_local_mpi_comm_world,ier)
#else
  recvbuf = sendbuf
#endif

  end subroutine max_all_all_i


!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_all_dp(sendbuf, recvbuf)

  use my_mpi_communicator
#ifdef USE_MPI
  use mpi
#endif

  implicit none

  double precision :: sendbuf, recvbuf

#ifdef USE_MPI
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MAX,my_local_mpi_comm_world,ier)
#else
  recvbuf = sendbuf
#endif

  end subroutine max_all_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine any_all_l(sendbuf, recvbuf)

  use my_mpi_communicator
#ifdef USE_MPI
  use mpi
#endif

  implicit none

  logical :: sendbuf, recvbuf

#ifdef USE_MPI
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_LOGICAL,MPI_LOR,my_local_mpi_comm_world,ier)
#else
  recvbuf = sendbuf
#endif

  end subroutine any_all_l

!
!-------------------------------------------------------------------------------------------------
!

! MPI summations

  subroutine sum_all_i(sendbuf, recvbuf)

  use my_mpi_communicator
#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer :: sendbuf, recvbuf

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,0,my_local_mpi_comm_world,ier)
#else
    recvbuf = sendbuf
#endif

  end subroutine sum_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_all_i(sendbuf, recvbuf)

  use my_mpi_communicator
#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer :: sendbuf, recvbuf

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,my_local_mpi_comm_world,ier)
#else
  recvbuf = sendbuf
#endif

  end subroutine sum_all_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_cr(sendbuf, recvbuf)

  use my_mpi_communicator
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

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_SUM,0,my_local_mpi_comm_world,ier)
#else
    recvbuf = sendbuf
#endif

  end subroutine sum_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_dp(sendbuf, recvbuf)

  use my_mpi_communicator
#ifdef USE_MPI
  use mpi
#endif

  implicit none

  double precision :: sendbuf, recvbuf

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,my_local_mpi_comm_world,ier)
#else
    recvbuf = sendbuf
#endif

  end subroutine sum_all_dp

!
!-------------------------------------------------------------------------------------------------
!


  subroutine sum_all_1Darray_i(sendbuf, recvbuf, nx)

  use my_mpi_communicator
#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer :: nx
  integer, dimension(nx) :: sendbuf, recvbuf

#ifdef USE_MPI
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,nx,MPI_INTEGER,MPI_SUM,0,my_local_mpi_comm_world,ier)
#else
  ! to avoid compiler warning
  recvbuf = sendbuf
  nx = nx
#endif

  end subroutine sum_all_1Darray_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_all_dp(sendbuf, recvbuf)

  use my_mpi_communicator
#ifdef USE_MPI
  use mpi
#endif

  implicit none

  double precision :: sendbuf, recvbuf

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_SUM,my_local_mpi_comm_world,ier)
#else
    recvbuf = sendbuf
#endif

  end subroutine sum_all_all_dp

!-------------------------------------------------------------------------------------------------
!
! Send/Receive MPI
!
!-------------------------------------------------------------------------------------------------

! asynchronuous send/receive

  subroutine isend_cr(sendbuf, sendcount, dest, sendtag, req)

  use my_mpi_communicator
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

  call MPI_ISEND(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag,my_local_mpi_comm_world,req,ier)
#else
  call stop_the_code('isend_cr not implemented for serial code')
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

 subroutine irecv_cr(recvbuf, recvcount, dest, recvtag, req)

 ! asynchronuous receive

  use my_mpi_communicator
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

  call MPI_IRECV(recvbuf,recvcount,CUSTOM_MPI_TYPE,dest,recvtag,my_local_mpi_comm_world,req,ier)
#else
  call stop_the_code('irecv_cr not implemented for serial code')
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


  subroutine send_singlei(sendbuf, dest, sendtag)

! synchronuous/blocking send

  use my_mpi_communicator
  use mpi

  implicit none

  integer :: dest,sendtag
  integer :: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf,1,MPI_INTEGER,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_singlei

#endif

!
!-------------------------------------------------------------------------------------------------
!
#ifdef USE_MPI

  subroutine send_singledp(sendbuf, dest, sendtag)

! synchronuous/blocking send

  use my_mpi_communicator
  use mpi

  implicit none

  integer :: dest,sendtag
  double precision :: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf,1,MPI_DOUBLE_PRECISION,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_singledp

#endif

!
!-------------------------------------------------------------------------------------------------
!
#ifdef USE_MPI

  subroutine recv_singlei(recvbuf, dest, recvtag )

! synchronuous/blocking receive
  use my_mpi_communicator
  use mpi

  implicit none

  integer :: dest,recvtag
  integer :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,1,MPI_INTEGER,dest,recvtag,my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_singlei

#endif

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_MPI

  subroutine recv_singledp(recvbuf, dest, recvtag )

! synchronuous/blocking receive
  use my_mpi_communicator
  use mpi

  implicit none

  integer :: dest,recvtag
  double precision :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,1,MPI_DOUBLE_PRECISION,dest,recvtag,my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_singledp

#endif

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_MPI


  subroutine send_i(sendbuf, sendcount, dest, sendtag)

! synchronuous/blocking send

  use my_mpi_communicator
  use mpi

  implicit none

  integer :: dest,sendtag
  integer :: sendcount
  integer,dimension(sendcount):: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf,sendcount,MPI_INTEGER,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_i

#endif

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_MPI

  subroutine send_dp(sendbuf, sendcount, dest, sendtag)

! synchronuous/blocking send

  use my_mpi_communicator
  use mpi

  implicit none

  integer :: dest,sendtag
  integer :: sendcount
  double precision,dimension(sendcount):: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf,sendcount,MPI_DOUBLE_PRECISION,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_dp

#endif

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_MPI

  subroutine recv_i(recvbuf,recvcount, dest, recvtag )

! synchronuous/blocking receive
  use my_mpi_communicator
  use mpi

  implicit none

  integer :: dest,recvtag
  integer :: recvcount
  integer,dimension(recvcount):: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,recvcount,MPI_INTEGER,dest,recvtag,my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_i

#endif

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_MPI

  subroutine recv_dp(recvbuf,recvcount, dest, recvtag )

! synchronuous/blocking receive
  use my_mpi_communicator
  use mpi

  implicit none

  integer :: dest,recvtag
  integer :: recvcount
  double precision,dimension(recvcount):: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,recvcount,MPI_DOUBLE_PRECISION,dest,recvtag,my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_dp

#endif


!-------------------------------------------------------------------------------------------------
!
! MPI gather helper
!
!-------------------------------------------------------------------------------------------------

  subroutine gather_all_i(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use my_mpi_communicator
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

  call MPI_GATHER(sendbuf,sendcnt,MPI_INTEGER,recvbuf,recvcount,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
#else
  recvbuf(:,0) = sendbuf(:)
#endif

  end subroutine gather_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_dp(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use my_mpi_communicator
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

  call MPI_GATHER(sendbuf,sendcnt,MPI_DOUBLE_PRECISION,recvbuf,recvcount,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_world,ier)
#else
  recvbuf(:,0) = sendbuf(:)
#endif

  end subroutine gather_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_singlei(sendbuf, recvbuf, NPROC)

  use my_mpi_communicator
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

  call MPI_GATHER(sendbuf,1,MPI_INTEGER,recvbuf,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
#else
  recvbuf(0) = sendbuf
#endif

  end subroutine gather_all_singlei

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_all_singlei(sendbuf, recvbuf, NPROC)

  use my_mpi_communicator
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

  call MPI_ALLGATHER(sendbuf,1,MPI_INTEGER,recvbuf,1,MPI_INTEGER,my_local_mpi_comm_world,ier)
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

  use my_mpi_communicator
#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer,intent(out) :: sizeval

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_COMM_SIZE(my_local_mpi_comm_world,sizeval,ier)
  if (ier /= 0 ) call stop_the_code('Error getting MPI world size')
#else
  ! single process
  sizeval = 1
#endif

  end subroutine world_size

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_rank(rank)

  use my_mpi_communicator
#ifdef USE_MPI
! standard include of the MPI library
  use mpi
#endif

  implicit none

  integer,intent(out) :: rank

#ifdef USE_MPI
  ! local parameters
  integer :: ier

  call MPI_COMM_RANK(my_local_mpi_comm_world,rank,ier)
  if (ier /= 0 ) call stop_the_code('Error getting MPI rank')
#else
  ! always returns master rank zero
  rank = 0
#endif

  end subroutine world_rank

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_duplicate(comm)

  use my_mpi_communicator
#ifdef USE_MPI
! standard include of the MPI library
  use mpi
#endif

  implicit none

  integer,intent(out) :: comm

#ifdef USE_MPI
  integer :: ier

  call MPI_COMM_DUP(my_local_mpi_comm_world,comm,ier)
  if (ier /= 0) call stop_the_code('error duplicating my_local_mpi_comm_world communicator')
#else
  comm = comm
#endif

  end subroutine world_duplicate

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_get_comm(comm)

  use my_mpi_communicator
#ifdef USE_MPI
! standard include of the MPI library
  use mpi
#endif

  implicit none

  integer,intent(out) :: comm

  comm = my_local_mpi_comm_world

  end subroutine world_get_comm

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_MPI
! create sub-communicators if needed, if running more than one earthquake from the same job.
! create a sub-communicator for each independent run;
! if there is a single run to do, then just copy the default communicator to the new one

  subroutine world_split()

  use my_mpi_communicator
! standard include of the MPI library
  use mpi

  use constants, only: MAX_STRING_LEN,OUTPUT_FILES, &
    IMAIN,ISTANDARD_OUTPUT,mygroup,I_should_read_the_database
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer :: sizeval,myrank,ier,key,my_group_for_bcast,my_local_rank_for_bcast,NPROC

  character(len=MAX_STRING_LEN) :: path_to_add

  if (NUMBER_OF_SIMULTANEOUS_RUNS <= 0) call stop_the_code('NUMBER_OF_SIMULTANEOUS_RUNS <= 0 makes no sense')

  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeval,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mod(sizeval,NUMBER_OF_SIMULTANEOUS_RUNS) /= 0) then
    if (myrank == 0) print *,'Error: the number of MPI processes ',sizeval, &
                            ' is not a multiple of NUMBER_OF_SIMULTANEOUS_RUNS = ',NUMBER_OF_SIMULTANEOUS_RUNS
    call stop_the_code( &
'the number of MPI processes is not a multiple of NUMBER_OF_SIMULTANEOUS_RUNS. Make sure you call meshfem2D with mpirun.')
  endif

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. IMAIN == ISTANDARD_OUTPUT) &
    call stop_the_code( &
'must not have IMAIN == ISTANDARD_OUTPUT when NUMBER_OF_SIMULTANEOUS_RUNS > 1 otherwise output to screen is mingled. &
                 & Change this in specfem/setup/constant.h.in and recompile.')

  if (NUMBER_OF_SIMULTANEOUS_RUNS == 1) then

    my_local_mpi_comm_world = MPI_COMM_WORLD

! no broadcast of the mesh and model databases to other runs in that case
    my_group_for_bcast = 0
    my_local_mpi_comm_for_bcast = MPI_COMM_NULL

  else

!--- create a subcommunicator for each independent run

    NPROC = sizeval / NUMBER_OF_SIMULTANEOUS_RUNS

!   create the different groups of processes, one for each independent run
    mygroup = myrank / NPROC
    key = myrank
    if (mygroup < 0 .or. mygroup > NUMBER_OF_SIMULTANEOUS_RUNS-1) call stop_the_code('invalid value of mygroup')

!   build the sub-communicators
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, mygroup, key, my_local_mpi_comm_world, ier)
    if (ier /= 0) call stop_the_code('error while trying to create the sub-communicators')

!   add the right directory for that run
!   (group numbers start at zero, but directory names start at run0001, thus we add one)
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    OUTPUT_FILES = path_to_add(1:len_trim(path_to_add))//OUTPUT_FILES(1:len_trim(OUTPUT_FILES))
    !IN_DATA_FILES = path_to_add(1:len_trim(path_to_add))//IN_DATA_FILES(1:len_trim(IN_DATA_FILES))

!--- create a subcommunicator to broadcast the identical mesh and model databases if needed
    if (BROADCAST_SAME_MESH_AND_MODEL) then

      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)
!     to broadcast the model, split along similar ranks per run instead
      my_group_for_bcast = mod(myrank,NPROC)
      key = myrank
      if (my_group_for_bcast < 0 .or. my_group_for_bcast > NPROC-1) call stop_the_code('invalid value of my_group_for_bcast')

!     build the sub-communicators
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, my_group_for_bcast, key, my_local_mpi_comm_for_bcast, ier)
      if (ier /= 0) call stop_the_code('error while trying to create the sub-communicators')

!     see if that process will need to read the mesh and model database and then broadcast it to others
      call MPI_COMM_RANK(my_local_mpi_comm_for_bcast,my_local_rank_for_bcast,ier)
      if (my_local_rank_for_bcast > 0) I_should_read_the_database = .false.

    else

! no broadcast of the mesh and model databases to other runs in that case
      my_group_for_bcast = 0
      my_local_mpi_comm_for_bcast = MPI_COMM_NULL

    endif

  endif

  end subroutine world_split

#endif
!
!-------------------------------------------------------------------------------------------------
!

! close sub-communicators if needed, if running more than one earthquake from the same job.
  subroutine world_unsplit()

  use my_mpi_communicator
#ifdef USE_MPI
! standard include of the MPI library
  use mpi
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer :: ier

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) then
    call MPI_COMM_FREE(my_local_mpi_comm_world,ier)
    if (BROADCAST_SAME_MESH_AND_MODEL) call MPI_COMM_FREE(my_local_mpi_comm_for_bcast,ier)
  endif
#endif

  end subroutine world_unsplit

