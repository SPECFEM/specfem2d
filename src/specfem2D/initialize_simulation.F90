
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
!               Pieyre Le Loher, pieyre DOT le-loher aT inria.fr
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


  subroutine initialize_simulation()

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par, only : nproc,myrank,ninterface_acoustic,ninterface_elastic,ninterface_poroelastic

  implicit none
  include "constants.h"

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
#else
  nproc = 1
  myrank = 0
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


!
!-------------------------------------------------------------------------------------------------
!


  subroutine initialize_simulation_domains()

  use specfem_par, only : any_acoustic,any_gravitoacoustic,any_elastic,any_poroelastic, &
                          anisotropic,acoustic,gravitoacoustic,elastic,poroelastic,porosity,anisotropy,kmato, &
                          nspec,nspec_allocate,p_sv,ATTENUATION_VISCOELASTIC_SOLID,count_nspec_acoustic
  implicit none
  include "constants.h"

  ! local parameters
  integer :: ispec

  ! initializes
  any_acoustic = .false.
  any_gravitoacoustic = .false.
  any_elastic = .false.
  any_poroelastic = .false.

  anisotropic(:) = .false.
  acoustic(:) = .false.
  gravitoacoustic(:) = .false.
  elastic(:) = .false.
  poroelastic(:) = .false.

  ! loops over all elements
  count_nspec_acoustic = 0
  do ispec = 1,nspec

    if( nint(porosity(kmato(ispec))) == 1 ) then
      ! assume acoustic domain
      ! if gravitoacoustic -> set by read_external_model
      acoustic(ispec) = .true.
      elastic(ispec) = .false.
      poroelastic(ispec) = .false.
      any_acoustic = .true.
      gravitoacoustic(ispec) = .false.
      any_gravitoacoustic = .false.
      count_nspec_acoustic = count_nspec_acoustic + 1
    else if( porosity(kmato(ispec)) < TINYVAL) then
      ! assume elastic domain
      elastic(ispec) = .true.
      poroelastic(ispec) = .false.
      any_elastic = .true.
      if(any(anisotropy(:,kmato(ispec)) /= 0)) then
         anisotropic(ispec) = .true.
      endif
    else
      ! assume poroelastic domain
      elastic(ispec) = .false.
      poroelastic(ispec) = .true.
      any_poroelastic = .true.
    endif

  enddo ! of do ispec = 1,nspec


  if(.not. p_sv .and. .not. any_elastic) then
    print *, '*************** WARNING ***************'
    print *, 'Surface (membrane) waves calculation needs an elastic medium'
    print *, '*************** WARNING ***************'
    stop
  endif
  if(.not. p_sv .and. (ATTENUATION_VISCOELASTIC_SOLID)) then
    print *, '*************** WARNING ***************'
    print *, 'Attenuation and anisotropy are not implemented for surface (membrane) waves calculation'
    print *, '*************** WARNING ***************'
    stop
  endif


  if(ATTENUATION_VISCOELASTIC_SOLID) then
    nspec_allocate = nspec
  else
    nspec_allocate = 1
  endif

  end subroutine initialize_simulation_domains
