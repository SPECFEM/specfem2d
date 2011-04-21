
!========================================================================
!
!                   S P E C F E M 2 D  Version 6.1
!                   ------------------------------
!
! Copyright Universite de Pau, CNRS and INRIA, France,
! and Princeton University / California Institute of Technology, USA.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
!               Pieyre Le Loher, pieyre DOT le-loher aT inria.fr
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
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


  subroutine initialize_simulation(nproc,myrank,NUMBER_OF_PASSES, &
                  ninterface_acoustic,ninterface_elastic,ninterface_poroelastic)

  implicit none
  include "constants.h"
#ifdef USE_MPI
  include "mpif.h"
#endif

  integer :: nproc,myrank,NUMBER_OF_PASSES
  integer :: ninterface_acoustic, ninterface_elastic,ninterface_poroelastic

  ! local parameters
  integer :: ier
  character(len=256)  :: prname

!***********************************************************************
!
!             i n i t i a l i z a t i o n    p h a s e
!
!***********************************************************************

#ifdef USE_MPI
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)
  if( ier /= 0 ) call exit_MPI('error MPI initialization')

  ! this is only used in the case of MPI because it distinguishes between inner and outer element
  ! in the MPI partitions, which is meaningless in the serial case
  if(FURTHER_REDUCE_CACHE_MISSES) then
    NUMBER_OF_PASSES = 2
  else
    NUMBER_OF_PASSES = 1
  endif

#else
  nproc = 1
  myrank = 0
  !ier = 0
  !ninterface_acoustic = 0
  !ninterface_elastic = 0
  !ninterface_poroelastic = 0
  !iproc = 0
  !ispec_inner = 0
  !ispec_outer = 0

  if(PERFORM_CUTHILL_MCKEE) then
    NUMBER_OF_PASSES = 2
  else
    NUMBER_OF_PASSES = 1
  endif
#endif

  ninterface_acoustic = 0
  ninterface_elastic = 0
  ninterface_poroelastic = 0

  ! determine if we write to file instead of standard output
  if(IOUT /= ISTANDARD_OUTPUT) then

#ifdef USE_MPI
    write(prname,240) myrank
 240 format('simulation_results',i5.5,'.txt')
#else
    prname = 'simulation_results.txt'
#endif

    open(IOUT,file=prname,status='unknown',action='write',iostat=ier)
    if( ier /= 0 ) call exit_MPI('error opening file simulation_results***.txt')

  endif

  end subroutine initialize_simulation
