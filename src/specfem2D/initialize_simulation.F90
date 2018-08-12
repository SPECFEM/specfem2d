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

  subroutine initialize_simulation()

#ifdef USE_MPI
  use mpi
#endif
  use constants, only: IMAIN,ISTANDARD_OUTPUT,SIZE_REAL,NSTAGE,OUTPUT_FILES
  use specfem_par
  use specfem_par_movie, only: cutsnaps

  implicit none

  include 'version.fh'

  ! local parameters
  integer :: ier

!***********************************************************************
!
!             i n i t i a l i z a t i o n    p h a s e
!
!***********************************************************************

  OUTPUT_FILES = trim(OUTPUT_FILES) ! TODO OUTPUT_FILES!

  ! number of MPI processes
  call world_size(NPROC)

  ! myrank is the rank of each process, between 0 and NPROC-1.
  ! as usual in MPI, process 0 is in charge of coordinating everything
  ! and also takes care of the main output
  call world_rank(myrank)

  ! check process setup
  if (NPROC < 1) call stop_the_code('should have NPROC >= 1')

  ! checks rank to make sure that myrank is zero for serial version
#ifdef USE_MPI
    if (myrank >= NPROC) call stop_the_code('Error: invalid MPI rank')
#else
    ! serial version: checks rank is initialized
    if (myrank /= 0) call stop_the_code('Error: process should have myrank zero for serial version')
#endif

  ! determine if we write to file instead of standard output
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) then
    ! sets main output file name
    ! opens for simulation output
    open(IMAIN,file=trim(OUTPUT_FILES)//'output_solver.txt',status='unknown',action='write',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(OUTPUT_FILES)//'output_solver.txt')
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)

#ifdef USE_MPI
    write(IMAIN,*) '**********************************************'
    write(IMAIN,*) '**** Specfem 2-D Solver - MPI version     ****'
    write(IMAIN,*) '**********************************************'
#else
    write(IMAIN,*) '**********************************************'
    write(IMAIN,*) '**** Specfem 2-D Solver - serial version  ****'
    write(IMAIN,*) '**********************************************'
#endif

    write(IMAIN,*)
    write(IMAIN,*) 'Running Git version of the code corresponding to ', git_commit_version
    write(IMAIN,*) 'dating ', git_date_version
    write(IMAIN,*)

#ifdef USE_MPI
    write(IMAIN,*) 'There are ',NPROC,' MPI processes'
    write(IMAIN,*) 'Processes are numbered from 0 to ',NPROC-1
    write(IMAIN,*)
    write(IMAIN,*) 'There is a total of ',NPROC,' slices'
    write(IMAIN,*)
#endif

    write(IMAIN,*)
    write(IMAIN,*) 'NDIM = ',NDIM
    write(IMAIN,*)
    write(IMAIN,*) 'NGLLX = ',NGLLX
    write(IMAIN,*) 'NGLLZ = ',NGLLZ
    write(IMAIN,*)
    ! write information about precision used for floating-point operations
    if (CUSTOM_REAL == SIZE_REAL) then
      write(IMAIN,*) 'using single precision for the calculations'
    else
      write(IMAIN,*) 'using double precision for the calculations'
    endif
    write(IMAIN,*)
    write(IMAIN,*) 'smallest and largest possible floating-point numbers are: ', &
                   tiny(1._CUSTOM_REAL),huge(1._CUSTOM_REAL)
    call flush_IMAIN()
  endif

  ! starts reading in Database file (header info, simulation flags, number of elements
  call read_mesh_for_init()

  ! checks flags
  call initialize_simulation_check()

  ! initializes GPU cards
  if (GPU_MODE) call initialize_GPU()

  ! ----------- initialization and determining parameters
  ! allocates mesh
  ! local to global indexing
  allocate(ibool(NGLLX,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating ibool array')

  ! mesh arrays
  allocate(xix(NGLLX,NGLLZ,nspec), &
           xiz(NGLLX,NGLLZ,nspec), &
           gammax(NGLLX,NGLLZ,nspec), &
           gammaz(NGLLX,NGLLZ,nspec), &
           jacobian(NGLLX,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating mesh arrays for databases')

  ! domain flags
  allocate(ispec_is_elastic(nspec), &
           ispec_is_acoustic(nspec), &
           ispec_is_poroelastic(nspec),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating domain flag arrays')

  ispec_is_elastic(:) = .false.
  ispec_is_acoustic(:) = .false.
  ispec_is_poroelastic(:) = .false.

  ! element property flags
  allocate(ispec_is_anisotropic(nspec), &
           ispec_is_PML(nspec), stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating element property flag arrays')

  ispec_is_anisotropic(:) = .false.
  ispec_is_PML(:) = .false.

  ! initializes domain interfaces
  ninterface_acoustic = 0
  ninterface_elastic = 0
  ninterface_poroelastic = 0

  ! time scheme
  ! defines number of stages of chosen time stepping scheme
  if (time_stepping_scheme == 1) then
    ! Newmark
    stage_time_scheme = 1
  else if (time_stepping_scheme == 2) then
    ! LDDRK
    stage_time_scheme = NSTAGE ! 6
  else if (time_stepping_scheme == 3) then
    ! Runge-Kutta
    stage_time_scheme = 4
  endif

  ! converts percentage
  cutsnaps = cutsnaps / 100.d0

  ! sets model flag
  if (trim(MODEL) == 'default') then
    assign_external_model = .false.
  else
    assign_external_model = .true.
  endif

  ! make sure NSTEP is a multiple of subsamp_seismos
  ! if not, increase it a little bit, to the next multiple
  if (mod(NSTEP,subsamp_seismos) /= 0) then
    NSTEP = (NSTEP/subsamp_seismos + 1)*subsamp_seismos

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'NSTEP is not a multiple of subsamp_seismos'
      write(IMAIN,*) 'thus increasing it automatically to the next multiple, which is ',NSTEP
      write(IMAIN,*)
    endif
  endif

  ! output seismograms at least once at the end of the simulation
  NSTEP_BETWEEN_OUTPUT_SEISMOS = min(NSTEP,NSTEP_BETWEEN_OUTPUT_SEISMOS)

  ! make sure NSTEP_BETWEEN_OUTPUT_SEISMOS is a multiple of subsamp_seismos
  if (mod(NSTEP_BETWEEN_OUTPUT_SEISMOS,subsamp_seismos) /= 0) then
    if (myrank == 0) then
      write(IMAIN,*) 'Invalid number of NSTEP_BETWEEN_OUTPUT_SEISMOS = ',NSTEP_BETWEEN_OUTPUT_SEISMOS
      write(IMAIN,*) 'Must be a multiple of subsamp_seismos = ',subsamp_seismos
    endif
    call stop_the_code('Error: NSTEP_BETWEEN_OUTPUT_SEISMOS must be a multiple of subsamp_seismos')
  endif

  end subroutine initialize_simulation


!
!-------------------------------------------------------------------------------------------------
!

  subroutine initialize_simulation_check()

  use constants, only: IMAIN,ISTANDARD_OUTPUT,USE_ENFORCE_FIELDS
  use specfem_par
  use specfem_par_movie

  implicit none

  ! synchronizes processes
  call synchronize_all()

  ! number of processes
  if (nproc_read_from_database < 1) call stop_the_code('should have nproc_read_from_database >= 1')

  ! check that the code is running with the requested nb of processes
  if (NPROC /= nproc_read_from_database) then
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Error: number of processors supposed to run on: ',nproc_read_from_database
      write(IMAIN,*) 'Error: number of MPI processors actually run on: ',NPROC
      write(IMAIN,*)
      if (IMAIN /= ISTANDARD_OUTPUT) then
        print *
        print *, 'Error specfem3D: number of processors supposed to run on: ',nproc_read_from_database
        print *, 'Error specfem3D: number of MPI processors actually run on: ',NPROC
        print *
      endif
    endif
    call exit_MPI(myrank,'wrong number of MPI processes, must always have NPROC == nproc_read_from_database')
  endif

  ! time scheme
  if (SIMULATION_TYPE == 3 .and. (time_stepping_scheme == 2 .or. time_stepping_scheme == 3)) &
    call stop_the_code('RK and LDDRK time scheme not supported for adjoint inversion')

  ! standard RK scheme
  if (time_stepping_scheme == 3) then
    if (NPROC > 1) &
      call stop_the_code('MPI support for standard Runge-Kutta scheme is not implemented yet')

    if (PML_BOUNDARY_CONDITIONS) &
      call stop_the_code('PML boundary conditions not implemented with standard Runge Kutta scheme yet')
  endif

  ! Bielak parameter setup
  if (add_Bielak_conditions .and. .not. initialfield) &
    call stop_the_code('need to have an initial field to add Bielak plane wave conditions')

  ! seismogram output
  if (seismotype < 1 .or. seismotype > 6) &
    call stop_the_code( &
'seismotype should be 1(=displ), 2(=veloc), 3(=accel), 4(=pressure), 5(=curl of displ) or 6(=the fluid potential)')

  if (SAVE_FORWARD .and. (seismotype /= 1 .and. seismotype /= 6)) then
    ! user warning
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '***** WARNING *****'
      write(IMAIN,*) 'SAVE_FORWARD simulation: uses seismotype = ',seismotype
      write(IMAIN,*) '  Note that seismograms usually must be in displacement/potential for (poro)elastic/acoustic domains'
      write(IMAIN,*) '  and seismotype should be 1 (elastic/poroelastic adjoint source) or 6 (acoustic adjoint source)'
      write(IMAIN,*) '*******************'
      write(IMAIN,*)
    endif
    ! safety stop
    !  stop
  endif

  ! image type
  if (imagetype_postscript < 1 .or. imagetype_postscript > 4) call stop_the_code('Wrong type for PostScript snapshots')

  ! checks attenuation setting
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. (ATTENUATION_PORO_FLUID_PART)) then
    print *, '*************** Error ***************'
    call stop_the_code('Anisotropy & Viscous damping are not presently implemented for adjoint calculations')
  endif

  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. &
      (ATTENUATION_VISCOELASTIC .or. ATTENUATION_VISCOACOUSTIC) &
      .and. (.not. UNDO_ATTENUATION_AND_OR_PML .and. .not. NO_BACKWARD_RECONSTRUCTION )) then
    print *, '*************** Error ***************'
    call stop_the_code(&
    'attenuation is only implemented for adjoint calculations with UNDO_ATTENUATION_AND_OR_PML or NO_BACKWARD_RECONSTRUCTION')
  endif

  if ((.not. P_SV) .and. ATTENUATION_VISCOELASTIC) then
    print *, '*************** WARNING ***************'
    print *, 'Attenuation and anisotropy are not implemented for surface (membrane) waves calculation'
    print *, '*************** WARNING ***************'
    call stop_the_code('Please set P_SV flag to .true. for simulations with attenuation and anisotropy')
  endif

  if (USE_ENFORCE_FIELDS .and. SIMULATION_TYPE /= 1) &
    call stop_the_code('USE_ENFORCE_FIELDS is not supported yet for adjoint/kernel simulations')

  ! synchronizes processes
  call synchronize_all()

  end subroutine initialize_simulation_check

!
!-------------------------------------------------------------------------------------------------
!

  subroutine initialize_GPU()

! initialization for GPU cards

  use constants, only: IMAIN
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
  if (NGLLX /= NGLLZ) &
    call stop_the_code('GPU mode can only be used if NGLLX == NGLLZ ')
  if (CUSTOM_REAL /= 4 ) &
    call stop_the_code('GPU mode runs only with CUSTOM_REAL == 4')
  if (PML_BOUNDARY_CONDITIONS) &
    call stop_the_code('for PML boundaries, GPU_MODE is not supported yet')

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

  ! synchronizes processes
  call synchronize_all()

  end subroutine initialize_GPU

