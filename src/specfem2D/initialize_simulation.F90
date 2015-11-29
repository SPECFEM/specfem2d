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

  use constants,only: IMAIN,ISTANDARD_OUTPUT

  use specfem_par, only : NPROC,myrank,GPU_MODE,ninterface_acoustic,ninterface_elastic,ninterface_poroelastic

  implicit none

  ! local parameters
  integer :: ier
  character(len=256)  :: prname

!***********************************************************************
!
!             i n i t i a l i z a t i o n    p h a s e
!
!***********************************************************************

  ! number of MPI processes
  call world_size(NPROC)

  ! myrank is the rank of each process, between 0 and NPROC-1.
  ! as usual in MPI, process 0 is in charge of coordinating everything
  ! and also takes care of the main output
  call world_rank(myrank)


  ! check process setup
  if (NPROC < 1) stop 'should have NPROC >= 1'

  ! determine if we write to file instead of standard output
  if (IMAIN /= ISTANDARD_OUTPUT) then
    ! sets main output file name
#ifdef USE_MPI
    write(prname,"('output_solver',i5.5,'.txt')") myrank
#else
    prname = 'output_solver.txt'
    ! serial version: checks rank is initialized
    if (myrank /= 0) stop 'process should have myrank zero'
#endif
    ! opens for simulation output
    open(IMAIN,file='OUTPUT_FILES/'//trim(prname),status='unknown',action='write',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file OUTPUT_FILES/output_solver***.txt')
  endif

  ! starts reading in Database file
  call read_databases_init()

  ! checks flags
  call initialize_simulation_check()

  ! initializes GPU cards
  if (GPU_MODE) call initialize_GPU()

  ! initializes domain interfaces
  ninterface_acoustic = 0
  ninterface_elastic = 0
  ninterface_poroelastic = 0

  end subroutine initialize_simulation


!
!-------------------------------------------------------------------------------------------------
!

  subroutine initialize_simulation_check()

  use specfem_par

  implicit none

  ! number of processes
  if (nproc_read_from_database < 1) stop 'should have nproc_read_from_database >= 1'

  ! checks if matching with mpi processes
  if (NPROC /= nproc_read_from_database) stop 'must always have nproc == nproc_read_from_database'

  ! time scheme
  if (SIMULATION_TYPE == 3 .and.(time_stepping_scheme == 2 .or. time_stepping_scheme == 3)) &
                                  stop 'RK and LDDRK time scheme not supported for adjoint inversion'

  end subroutine initialize_simulation_check

!
!-------------------------------------------------------------------------------------------------
!


  subroutine initialize_simulation_domains()

  use constants,only: TINYVAL

  use specfem_par, only : any_acoustic,any_gravitoacoustic,any_elastic,any_poroelastic, &
    ispec_is_anisotropic,ispec_is_acoustic,ispec_is_gravitoacoustic,ispec_is_elastic,ispec_is_poroelastic, &
    porosity,anisotropy,kmato, &
    nspec,nspec_allocate,p_sv,ATTENUATION_VISCOELASTIC_SOLID,count_nspec_acoustic

  implicit none

  ! local parameters
  integer :: ispec

  ! initializes
  any_acoustic = .false.
  any_gravitoacoustic = .false.
  any_elastic = .false.
  any_poroelastic = .false.

  ispec_is_anisotropic(:) = .false.
  ispec_is_acoustic(:) = .false.
  ispec_is_gravitoacoustic(:) = .false.
  ispec_is_elastic(:) = .false.
  ispec_is_poroelastic(:) = .false.

  ! loops over all elements
  count_nspec_acoustic = 0
  do ispec = 1,nspec

    ! checks domain properties
    if (nint(porosity(kmato(ispec))) == 1) then
      ! assume acoustic domain
      ! if gravitoacoustic -> set by read_external_model
      ispec_is_acoustic(ispec) = .true.
      ispec_is_elastic(ispec) = .false.
      ispec_is_poroelastic(ispec) = .false.
      any_acoustic = .true.
      ispec_is_gravitoacoustic(ispec) = .false.
      any_gravitoacoustic = .false.
      count_nspec_acoustic = count_nspec_acoustic + 1

    else if (porosity(kmato(ispec)) < TINYVAL) then
      ! assume elastic domain
      ispec_is_elastic(ispec) = .true.
      ispec_is_poroelastic(ispec) = .false.
      any_elastic = .true.
      if (any(anisotropy(:,kmato(ispec)) /= 0)) then
        ispec_is_anisotropic(ispec) = .true.
      endif

    else
      ! assume poroelastic domain
      ispec_is_elastic(ispec) = .false.
      ispec_is_poroelastic(ispec) = .true.
      any_poroelastic = .true.
    endif

  enddo ! of do ispec = 1,nspec


  if (.not. p_sv .and. .not. any_elastic) then
    print *, '*************** WARNING ***************'
    print *, 'Surface (membrane) waves calculation needs an elastic medium'
    print *, '*************** WARNING ***************'
    stop
  endif
  if (.not. p_sv .and. (ATTENUATION_VISCOELASTIC_SOLID)) then
    print *, '*************** WARNING ***************'
    print *, 'Attenuation and anisotropy are not implemented for surface (membrane) waves calculation'
    print *, '*************** WARNING ***************'
    stop
  endif


  if (ATTENUATION_VISCOELASTIC_SOLID) then
    nspec_allocate = nspec
  else
    nspec_allocate = 1
  endif

  end subroutine initialize_simulation_domains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine initialize_GPU()

! initialization for GPU cards

  use specfem_par

  implicit none

  ! local parameters
  integer :: ncuda_devices,ncuda_devices_min,ncuda_devices_max

  ! GPU_MODE now defined in Par_file
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) "GPU_MODE Active."
    call flush_IMAIN()
  endif

  ! check for GPU runs
  if (NGLLX /= 5 .or. NGLLZ /= 5 ) &
    stop 'GPU mode can only be used if NGLLX == NGLLZ == 5'
  if (CUSTOM_REAL /= 4 ) &
    stop 'GPU mode runs only with CUSTOM_REAL == 4'

  ! initializes GPU and outputs info to files for all processes
  call initialize_cuda_device(myrank,ncuda_devices)

  ! synchronizes all processes
  call synchronize_all()

  ! collects min/max of local devices found for statistics
#ifdef USE_MPI
  call min_all_i(ncuda_devices,ncuda_devices_min)
  call max_all_i(ncuda_devices,ncuda_devices_max)
#else
  ncuda_devices_min = ncuda_devices
  ncuda_devices_max = ncuda_devices
#endif

  if (myrank == 0) then
    write(IMAIN,*) "GPU number of devices per node: min =",ncuda_devices_min
    write(IMAIN,*) "                                max =",ncuda_devices_max
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine initialize_GPU

