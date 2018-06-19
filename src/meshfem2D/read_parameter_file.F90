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

  subroutine read_parameter_file(imesher,BROADCAST_AFTER_READ)

! reads in DATA/Par_file

  use constants, only: IMAIN
  use shared_parameters

  implicit none

  integer, intent(in) :: imesher
  logical, intent(in) :: BROADCAST_AFTER_READ

  ! initializes
  call read_parameter_file_init()

  ! only master process reads in Par_file
  if (myrank == 0) then
    ! opens file Par_file
    call open_parameter_file()

    ! reads only parameters (without receiver-line section, material tables or region definitions)
    call read_parameter_file_only()

    ! reads receiver lines
    call read_parameter_file_receiversets()

    ! reads material definitions
    call read_material_table()

    ! mesher reads in internal region table for setting up mesh elements
    if (imesher == 1 .and. (.not. read_external_mesh) ) then
      ! internal meshing
      ! user output
      write(IMAIN,*)
      write(IMAIN,*) 'Mesh from internal meshing:'
      ! reads interface definitions from interface file (we need to have nxread & nzread value for checking regions)
      call read_interfaces_file()

      ! internal meshing
      nx = nxread
      nz = nzread

      ! setup mesh array
      ! multiply by 2 if elements have 9 nodes
      if (ngnod == 9) then
        nx = nx * 2
        nz = nz * 2
        nz_layer(:) = nz_layer(:) * 2
      endif

      ! total number of elements
      nelmnts = nxread * nzread

      ! reads material regions defined in Par_file
      call read_regions()
    endif

    ! closes file Par_file
    call close_parameter_file()
  endif

  ! master process broadcasts to all
  ! note: this is only needed at the moment for the solver to setup a simulation run
  if (BROADCAST_AFTER_READ) then
    call bcast_all_singlei(SIMULATION_TYPE)
    call bcast_all_singlei(NOISE_TOMOGRAPHY)
    call bcast_all_singlel(SAVE_FORWARD)

    call bcast_all_singlei(NPROC)
    call bcast_all_singlei(partitioning_method)
    call bcast_all_singlei(ngnod)

    call bcast_all_singlei(NSTEP)
    call bcast_all_singledp(DT)

    call bcast_all_singlei(time_stepping_scheme)
    call bcast_all_singlel(AXISYM)
    call bcast_all_singlel(P_SV)
    call bcast_all_singlel(GPU_MODE)
    call bcast_all_singlei(setup_with_binary_database)

    call bcast_all_string(MODEL)
    call bcast_all_string(SAVE_MODEL)

    call bcast_all_singlel(ATTENUATION_VISCOELASTIC)
    call bcast_all_singlel(ATTENUATION_VISCOACOUSTIC)
    call bcast_all_singlei(N_SLS)
    call bcast_all_singledp(ATTENUATION_f0_REFERENCE)
    call bcast_all_singlel(READ_VELOCITIES_AT_f0)
    call bcast_all_singlel(USE_SOLVOPT)
    call bcast_all_singlel(ATTENUATION_PORO_FLUID_PART)
    call bcast_all_singledp(Q0_poroelastic)
    call bcast_all_singledp(freq0_poroelastic)

    call bcast_all_singlel(UNDO_ATTENUATION_AND_OR_PML)
    call bcast_all_singlei(NT_DUMP_ATTENUATION)

    call bcast_all_singlei(NSOURCES)
    call bcast_all_singlel(force_normal_to_surface)
    call bcast_all_singlel(initialfield)
    call bcast_all_singlel(add_Bielak_conditions_bottom)
    call bcast_all_singlel(add_Bielak_conditions_right)
    call bcast_all_singlel(add_Bielak_conditions_top)
    call bcast_all_singlel(add_Bielak_conditions_left)
    call bcast_all_singlel(ACOUSTIC_FORCING)

    call bcast_all_singlei(seismotype)
    call bcast_all_singlei(subsamp_seismos)
    call bcast_all_singlel(USE_TRICK_FOR_BETTER_PRESSURE)
    call bcast_all_singlei(NSTEP_BETWEEN_OUTPUT_SEISMOS)
    call bcast_all_singlel(COMPUTE_INTEGRATED_ENERGY_FIELD)
    call bcast_all_singledp(USER_T0)
    call bcast_all_singlel(save_ASCII_seismograms)
    call bcast_all_singlel(save_binary_seismograms_single)
    call bcast_all_singlel(save_binary_seismograms_double)
    call bcast_all_singlel(SU_FORMAT)
    call bcast_all_singlel(use_existing_STATIONS)
    call bcast_all_singlei(nreceiversets)
    call bcast_all_singledp(anglerec)
    call bcast_all_singlel(rec_normal_to_surface)

    call bcast_all_singlel(save_ASCII_kernels)
    call bcast_all_singlei(NSTEP_BETWEEN_COMPUTE_KERNELS)
    call bcast_all_singlel(NO_BACKWARD_RECONSTRUCTION)


    call bcast_all_singlel(STACEY_ABSORBING_CONDITIONS)
    call bcast_all_singlel(ADD_PERIODIC_CONDITIONS)
    call bcast_all_singledp(PERIODIC_HORIZ_DIST)

    call bcast_all_singlei(nbmodels)
    call bcast_all_string(TOMOGRAPHY_FILE)
    call bcast_all_singlel(read_external_mesh)

    if (.not. read_external_mesh) then
      call bcast_all_singlel(absorbbottom)
      call bcast_all_singlel(absorbright)
      call bcast_all_singlel(absorbtop)
      call bcast_all_singlel(absorbleft)
    endif

    call bcast_all_singlei(NSTEP_BETWEEN_OUTPUT_INFO)
    call bcast_all_singlel(output_grid_Gnuplot)
    call bcast_all_singlel(output_grid_ASCII)
    call bcast_all_singlel(OUTPUT_ENERGY)
    call bcast_all_singlei(NTSTEP_BETWEEN_OUTPUT_ENERGY)

    call bcast_all_singlei(NSTEP_BETWEEN_OUTPUT_IMAGES)
    call bcast_all_singledp(cutsnaps)

    call bcast_all_singlel(output_color_image)
    call bcast_all_singlei(imagetype_JPEG)
    call bcast_all_singledp(factor_subsample_image)
    call bcast_all_singlel(USE_CONSTANT_MAX_AMPLITUDE)
    call bcast_all_singledp(CONSTANT_MAX_AMPLITUDE_TO_USE)
    call bcast_all_singledp(POWER_DISPLAY_COLOR)
    call bcast_all_singlel(DRAW_SOURCES_AND_RECEIVERS)
    call bcast_all_singlel(DRAW_WATER_IN_BLUE)
    call bcast_all_singlel(USE_SNAPSHOT_NUMBER_IN_FILENAME)

    call bcast_all_singlel(output_postscript_snapshot)
    call bcast_all_singlei(imagetype_postscript)
    call bcast_all_singlel(meshvect)
    call bcast_all_singlel(modelvect)
    call bcast_all_singlel(boundvect)
    call bcast_all_singlel(interpol)
    call bcast_all_singlei(pointsdisp)
    call bcast_all_singlei(subsamp_postscript)
    call bcast_all_singledp(sizemax_arrows)
    call bcast_all_singlel(US_LETTER)

    call bcast_all_singlel(output_wavefield_dumps)
    call bcast_all_singlei(imagetype_wavefield_dumps)
    call bcast_all_singlel(use_binary_for_wavefield_dumps)

    call bcast_all_singlei(NUMBER_OF_SIMULTANEOUS_RUNS)
    call bcast_all_singlel(BROADCAST_SAME_MESH_AND_MODEL)
  endif

  ! derive additional settings/flags based on input parameters
  call read_parameter_file_derive_flags()

  end subroutine read_parameter_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_parameter_file_init()

! initializes the variables

  use shared_parameters

  implicit none

  ! external meshing
  mesh_file = ''
  nodes_coords_file = ''
  materials_file = ''
  free_surface_file = ''
  axial_elements_file = ''
  absorbing_surface_file = ''
  acoustic_forcing_surface_file = ''
  absorbing_cpml_file = ''
  tangential_detection_curve_file = ''

  ! internal meshing
  interfacesfile = ''
  xmin_param = 0.d0
  xmax_param = 0.d0
  nx_param = 0

  absorbbottom = .false.
  absorbright = .false.
  absorbtop = .false.
  absorbleft = .false.

  nbregions = 0

  end subroutine read_parameter_file_init

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_parameter_file_only()

! reads only parameters without receiver-line section and material tables

  use constants, only: IMAIN
  use shared_parameters

  implicit none

  ! local parameters
  integer :: i,irange
  logical :: some_parameters_missing_from_Par_file

  integer, external :: err_occurred

!! DK DK to detect discontinued parameters
  double precision :: f0_attenuation

  !--------------------------------------------------------------------
  !
  ! simulation input paramters
  !
  !--------------------------------------------------------------------

  some_parameters_missing_from_Par_file = .false.

  ! read file names and path for output
  call read_value_string_p(title, 'title')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'title                           = Title of my simulation'
    write(*,*)
  endif

  ! read type of simulation
  call read_value_integer_p(SIMULATION_TYPE, 'SIMULATION_TYPE')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'SIMULATION_TYPE                 = 1'
    write(*,*)
  endif

  call read_value_integer_p(NOISE_TOMOGRAPHY, 'NOISE_TOMOGRAPHY')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'NOISE_TOMOGRAPHY                = 0'
    write(*,*)
  endif

  call read_value_logical_p(SAVE_FORWARD, 'SAVE_FORWARD')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'SAVE_FORWARD                    = .false.'
    write(*,*)
  endif

  ! read info about partitioning
  call read_value_integer_p(NPROC, 'NPROC')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'NPROC                           = 1'
    write(*,*)
  endif

  call read_value_integer_p(partitioning_method, 'partitioning_method')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'partitioning_method             = 3'
    write(*,*)
  endif

  call read_value_integer_p(ngnod, 'ngnod')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'ngnod                           = 9'
    write(*,*)
  endif

  ! read time step parameters
  call read_value_integer_p(NSTEP, 'NSTEP')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'NSTEP                           = 3000'
    write(*,*)
  endif

  call read_value_double_precision_p(DT, 'DT')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'DT                              = 1.0e-3'
    write(*,*)
  endif

  call read_value_integer_p(time_stepping_scheme, 'time_stepping_scheme')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'time_stepping_scheme            = 1'
    write(*,*)
  endif

  ! axisymmetric (2.5D) or Cartesian planar (2D) simulation
  call read_value_logical_p(AXISYM, 'AXISYM')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'AXISYM                          = .false.'
    write(*,*)
  endif

  ! determine if body or surface (membrane) waves calculation
  call read_value_logical_p(P_SV, 'P_SV')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'P_SV                            = .true.'
    write(*,*)
  endif

  call read_value_logical_p(GPU_MODE, 'GPU_MODE')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'GPU_MODE                        = .false.'
    write(*,*)
  endif

  call read_value_integer_p(setup_with_binary_database, 'setup_with_binary_database')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'setup_with_binary_database      = 0'
    write(*,*)
  endif

  call read_value_string_p(MODEL, 'MODEL')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'MODEL                           = default'
    write(*,*)
  endif

  call read_value_string_p(SAVE_MODEL, 'SAVE_MODEL')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'SAVE_MODEL                      = default'
    write(*,*)
  endif

  !--------------------------------------------------------------------
  !
  ! attenuation
  !
  !--------------------------------------------------------------------

  call read_value_logical_p(ATTENUATION_VISCOELASTIC, 'ATTENUATION_VISCOELASTIC')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'ATTENUATION_VISCOELASTIC        = .false.'
    write(*,*)
  endif

  ! read viscous attenuation parameters (acoustic media)
  call read_value_logical_p(ATTENUATION_VISCOACOUSTIC, 'ATTENUATION_VISCOACOUSTIC')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'ATTENUATION_VISCOACOUSTIC       = .false.'
    write(*,*)
  endif

  ! read constants for attenuation
  call read_value_integer_p(N_SLS, 'N_SLS')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'N_SLS                           = 3'
    write(*,*)
  endif

  call read_value_double_precision_p(ATTENUATION_f0_REFERENCE, 'ATTENUATION_f0_REFERENCE')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'ATTENUATION_f0_REFERENCE        = 5.196'
    write(*,*)
  endif

!! DK DK discontinued parameter
  call read_value_double_precision_p(f0_attenuation, 'f0_attenuation')
! if this parameter exists in the Par_file
  if (err_occurred() == 0) then
    write(*,'(a)') 'Parameter f0_attenuation in the Par_file is now called ATTENUATION_f0_REFERENCE'
    write(*,'(a)') 'in order to use the same name as in the 3D code (SPECFEM3D).'
    write(*,'(a)') 'Please rename it in your Par_file and start the code again.'
    write(*,*)
    call stop_the_code('Error: parameter f0_attenuation should be renamed in your Par_file')
  endif

  call read_value_logical_p(READ_VELOCITIES_AT_f0, 'READ_VELOCITIES_AT_f0')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'READ_VELOCITIES_AT_f0           = .false.'
    write(*,*)
  endif

  call read_value_logical_p(USE_SOLVOPT, 'USE_SOLVOPT')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'USE_SOLVOPT                     = .false.'
    write(*,*)
  endif

  ! read viscous attenuation parameters (poroelastic media)
  call read_value_logical_p(ATTENUATION_PORO_FLUID_PART, 'ATTENUATION_PORO_FLUID_PART')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'ATTENUATION_PORO_FLUID_PART     = .false.'
    write(*,*)
  endif

  call read_value_double_precision_p(Q0_poroelastic, 'Q0_poroelastic')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'Q0_poroelastic                  = 1'
    write(*,*)
  endif

  call read_value_double_precision_p(freq0_poroelastic, 'freq0_poroelastic')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'freq0_poroelastic               = 10'
    write(*,*)
  endif

  call read_value_logical_p(UNDO_ATTENUATION_AND_OR_PML, 'UNDO_ATTENUATION_AND_OR_PML')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'UNDO_ATTENUATION_AND_OR_PML     = .false.'
    write(*,*)
  endif

  call read_value_integer_p(NT_DUMP_ATTENUATION, 'NT_DUMP_ATTENUATION')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'NT_DUMP_ATTENUATION             = 500'
    write(*,*)
  endif

  !--------------------------------------------------------------------
  !
  ! sources
  !
  !--------------------------------------------------------------------

  ! read source infos
  call read_value_integer_p(NSOURCES, 'NSOURCES')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'NSOURCES                        = 1'
    write(*,*)
  endif

  call read_value_logical_p(force_normal_to_surface, 'force_normal_to_surface')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'force_normal_to_surface         = .false.'
    write(*,*)
  endif

  call read_value_logical_p(initialfield, 'initialfield')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'initialfield                    = .false.'
    write(*,*)
  endif

  call read_value_logical_p(add_Bielak_conditions_bottom, 'add_Bielak_conditions_bottom')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'add_Bielak_conditions_bottom    = .false.'
    write(*,*)
  endif

  call read_value_logical_p(add_Bielak_conditions_right, 'add_Bielak_conditions_right')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'add_Bielak_conditions_right     = .false.'
    write(*,*)
  endif

  call read_value_logical_p(add_Bielak_conditions_top, 'add_Bielak_conditions_top')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'add_Bielak_conditions_top       = .false.'
    write(*,*)
  endif

  call read_value_logical_p(add_Bielak_conditions_left, 'add_Bielak_conditions_left')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'add_Bielak_conditions_left      = .false.'
    write(*,*)
  endif

  ! read acoustic forcing flag
  call read_value_logical_p(ACOUSTIC_FORCING, 'ACOUSTIC_FORCING')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'ACOUSTIC_FORCING                = .false.'
    write(*,*)
  endif

  !--------------------------------------------------------------------
  !
  ! receivers
  !
  !--------------------------------------------------------------------

  ! read receiver line parameters
  call read_value_integer_p(seismotype, 'seismotype')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'seismotype                      = 1'
    write(*,*)
  endif

  call read_value_integer_p(subsamp_seismos, 'subsamp_seismos')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'subsamp_seismos                 = 1'
    write(*,*)
  endif

  call read_value_logical_p(USE_TRICK_FOR_BETTER_PRESSURE, 'USE_TRICK_FOR_BETTER_PRESSURE')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'USE_TRICK_FOR_BETTER_PRESSURE   = .false.'
    write(*,*)
  endif

  call read_value_integer_p(NSTEP_BETWEEN_OUTPUT_SEISMOS, 'NSTEP_BETWEEN_OUTPUT_SEISMOS')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'NSTEP_BETWEEN_OUTPUT_SEISMOS    = 1000'
    write(*,*)
  endif

  call read_value_double_precision_p(USER_T0, 'USER_T0')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'USER_T0                         = 0.0d0'
    write(*,*)
  endif

  call read_value_logical_p(save_ASCII_seismograms, 'save_ASCII_seismograms')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'save_ASCII_seismograms          = .true.'
    write(*,*)
  endif

  call read_value_logical_p(save_binary_seismograms_single, 'save_binary_seismograms_single')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'save_binary_seismograms_single  = .true.'
    write(*,*)
  endif

  call read_value_logical_p(save_binary_seismograms_double, 'save_binary_seismograms_double')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'save_binary_seismograms_double  = .false.'
    write(*,*)
  endif

  call read_value_logical_p(SU_FORMAT, 'SU_FORMAT')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'SU_FORMAT                       = .false.'
    write(*,*)
  endif

  call read_value_logical_p(use_existing_STATIONS, 'use_existing_STATIONS')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'use_existing_STATIONS           = .false.'
    write(*,*)
  endif

  call read_value_integer_p(nreceiversets, 'nreceiversets')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'nreceiversets                   = 1'
    write(*,*)
  endif

  call read_value_double_precision_p(anglerec, 'anglerec')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'anglerec                        = 0.d0'
    write(*,*)
  endif

  call read_value_logical_p(rec_normal_to_surface, 'rec_normal_to_surface')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'rec_normal_to_surface           = .false.'
    write(*,*)
  endif

  ! receiver sets will be read in later...

  !--------------------------------------------------------------------
  !
  ! adjoint kernel
  !
  !--------------------------------------------------------------------

  call read_value_logical_p(save_ASCII_kernels, 'save_ASCII_kernels')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'save_ASCII_kernels              = .true.'
    write(*,*)
  endif

  call read_value_integer_p(NSTEP_BETWEEN_COMPUTE_KERNELS, 'NSTEP_BETWEEN_COMPUTE_KERNELS')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'NSTEP_BETWEEN_COMPUTE_KERNELS             = 1'
    write(*,*)
  endif

  call read_value_logical_p(NO_BACKWARD_RECONSTRUCTION,'NO_BACKWARD_RECONSTRUCTION')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'NO_BACKWARD_RECONSTRUCTIO           = .false.'
    write(*,*)
  endif


  !--------------------------------------------------------------------
  !
  ! boundary conditions
  !
  !--------------------------------------------------------------------

  call read_value_logical_p(PML_BOUNDARY_CONDITIONS, 'PML_BOUNDARY_CONDITIONS')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'PML_BOUNDARY_CONDITIONS         = .true.'
    write(*,*)
  endif

  call read_value_integer_p(NELEM_PML_THICKNESS, 'NELEM_PML_THICKNESS')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'NELEM_PML_THICKNESS             = 3'
    write(*,*)
  endif

  call read_value_logical_p(ROTATE_PML_ACTIVATE, 'ROTATE_PML_ACTIVATE')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'ROTATE_PML_ACTIVATE             = .false.'
    write(*,*)
  endif

  call read_value_double_precision_p(ROTATE_PML_ANGLE, 'ROTATE_PML_ANGLE')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'ROTATE_PML_ANGLE                = 30.'
    write(*,*)
  endif

  call read_value_double_precision_p(K_MIN_PML, 'K_MIN_PML')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'K_MIN_PML                       = 1.d0'
    write(*,*)
  endif

  call read_value_double_precision_p(K_MAX_PML, 'K_MAX_PML')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'K_MAX_PML                       = 1.d0'
    write(*,*)
  endif

  call read_value_double_precision_p(damping_change_factor_acoustic, 'damping_change_factor_acoustic')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'damping_change_factor_acoustic  = 0.5d0'
    write(*,*)
  endif

  call read_value_double_precision_p(damping_change_factor_elastic, 'damping_change_factor_elastic')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'damping_change_factor_elastic   = 1.d0'
    write(*,*)
  endif

  call read_value_logical_p(PML_PARAMETER_ADJUSTMENT, 'PML_PARAMETER_ADJUSTMENT')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'PML_PARAMETER_ADJUSTMENT        = .false.'
    write(*,*)
  endif

  ! boolean defining whether to use any absorbing boundaries
  call read_value_logical_p(STACEY_ABSORBING_CONDITIONS, 'STACEY_ABSORBING_CONDITIONS')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'STACEY_ABSORBING_CONDITIONS     = .false.'
    write(*,*)
  endif

  call read_value_logical_p(ADD_PERIODIC_CONDITIONS, 'ADD_PERIODIC_CONDITIONS')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'ADD_PERIODIC_CONDITIONS         = .false.'
    write(*,*)
  endif

  call read_value_double_precision_p(PERIODIC_HORIZ_DIST, 'PERIODIC_HORIZ_DIST')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'PERIODIC_HORIZ_DIST             = 4000.'
    write(*,*)
  endif

  !--------------------------------------------------------------------
  !
  ! velocity and density models
  !
  !--------------------------------------------------------------------

  ! read the different material materials (i.e. the number of models)
  call read_value_integer_p(nbmodels, 'nbmodels')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'nbmodels                        = 1'
    write(*,*)
  endif

  ! material definitions will be read later on...

  call read_value_string_p(TOMOGRAPHY_FILE, 'TOMOGRAPHY_FILE')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'TOMOGRAPHY_FILE                 = ./DATA/tomo_file.xyz'
    write(*,*)
  endif

  ! boolean defining whether internal or external mesh
  call read_value_logical_p(read_external_mesh, 'read_external_mesh')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'read_external_mesh              = .false.'
    write(*,*)
  endif

  !--------------------------------------------------------------------
  !
  ! parameters external / internal meshing
  !
  !--------------------------------------------------------------------

  !-----------------
  ! external mesh parameters

  if (read_external_mesh) then

    ! read info about external mesh
    call read_value_string_p(mesh_file, 'mesh_file')
    if (err_occurred() /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'mesh_file                       = ./DATA/mesh_file'
      write(*,*)
    endif

    call read_value_string_p(nodes_coords_file, 'nodes_coords_file')
    if (err_occurred() /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'nodes_coords_file               = ./DATA/nodes_coords_file'
      write(*,*)
    endif

    call read_value_string_p(materials_file, 'materials_file')
    if (err_occurred() /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'materials_file                  = ./DATA/materials_file'
      write(*,*)
    endif

    call read_value_string_p(free_surface_file, 'free_surface_file')
    if (err_occurred() /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'free_surface_file               = ./DATA/free_surface_file'
      write(*,*)
    endif

    call read_value_string_p(axial_elements_file, 'axial_elements_file')
    if (err_occurred() /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'axial_elements_file             = ./DATA/axial_elements_file'
      write(*,*)
    endif

    call read_value_string_p(absorbing_surface_file, 'absorbing_surface_file')
    if (err_occurred() /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'absorbing_surface_file          = ./DATA/absorbing_surface_file'
      write(*,*)
    endif

    call read_value_string_p(acoustic_forcing_surface_file, 'acoustic_forcing_surface_file')
    if (err_occurred() /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'acoustic_forcing_surface_file   = ./DATA/MSH/Surf_acforcing_Bottom_enforcing_mesh'
      write(*,*)
    endif

    call read_value_string_p(absorbing_cpml_file, 'absorbing_cpml_file')
    if (err_occurred() /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'absorbing_cpml_file             = ./DATA/absorbing_cpml_file'
      write(*,*)
    endif

    call read_value_string_p(tangential_detection_curve_file, 'tangential_detection_curve_file')
    if (err_occurred() /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'tangential_detection_curve_file = ./DATA/courbe_eros_nodes'
      write(*,*)
    endif

  else

    !-----------------
    ! internal mesh parameters

    ! interfaces file
    call read_value_string_p(interfacesfile, 'interfacesfile')
    if (err_occurred() /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'interfacesfile                  = DATA/interfaces_simple_topo_curved.dat'
      write(*,*)
    endif

    ! read grid parameters
    call read_value_double_precision_p(xmin_param, 'xmin')
    if (err_occurred() /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'xmin                            = 0.d0'
      write(*,*)
    endif

    call read_value_double_precision_p(xmax_param, 'xmax')
    if (err_occurred() /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'xmax                            = 4000.d0'
      write(*,*)
    endif

    call read_value_integer_p(nx_param, 'nx')
    if (err_occurred() /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'nx                              = 80'
      write(*,*)
    endif

    ! read absorbing boundary parameters
    call read_value_logical_p(absorbbottom, 'absorbbottom')
    if (err_occurred() /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'absorbbottom                    = .true.'
      write(*,*)
    endif

    call read_value_logical_p(absorbright, 'absorbright')
    if (err_occurred() /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'absorbright                     = .true.'
      write(*,*)
    endif

    call read_value_logical_p(absorbtop, 'absorbtop')
    if (err_occurred() /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'absorbtop                       = .false.'
      write(*,*)
    endif

    call read_value_logical_p(absorbleft, 'absorbleft')
    if (err_occurred() /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'absorbleft                      = .true.'
      write(*,*)
    endif

    call read_value_integer_p(nbregions, 'nbregions')
    if (err_occurred() /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'nbregions                       = 1'
      write(*,*)
    endif
! note: if internal mesh, then region tables will be read in by read_regions (from meshfem2D)
  endif


  !--------------------------------------------------------------------
  !
  ! display parameters
  !
  !--------------------------------------------------------------------

  ! read display parameters
  call read_value_integer_p(NSTEP_BETWEEN_OUTPUT_INFO, 'NSTEP_BETWEEN_OUTPUT_INFO')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'NSTEP_BETWEEN_OUTPUT_INFO       = 100'
    write(*,*)
  endif

  call read_value_logical_p(output_grid_Gnuplot, 'output_grid_Gnuplot')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'output_grid_Gnuplot             = .false.'
    write(*,*)
  endif

  call read_value_logical_p(output_grid_ASCII, 'output_grid_ASCII')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'output_grid_ASCII               = .false.'
    write(*,*)
  endif

  call read_value_logical_p(OUTPUT_ENERGY, 'OUTPUT_ENERGY')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'OUTPUT_ENERGY                   = .false.'
    write(*,*)
  endif

  call read_value_integer_p(NTSTEP_BETWEEN_OUTPUT_ENERGY, 'NTSTEP_BETWEEN_OUTPUT_ENERGY')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'NTSTEP_BETWEEN_OUTPUT_ENERGY    = 10'
    write(*,*)
  endif

  call read_value_logical_p(COMPUTE_INTEGRATED_ENERGY_FIELD, 'COMPUTE_INTEGRATED_ENERGY_FIELD')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'COMPUTE_INTEGRATED_ENERGY_FIELD = .false.'
    write(*,*)
  endif

  !--------------------------------------------------------------------
  !
  ! movies/images/snapshots
  !
  !--------------------------------------------------------------------

  call read_value_integer_p(NSTEP_BETWEEN_OUTPUT_IMAGES, 'NSTEP_BETWEEN_OUTPUT_IMAGES')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'NSTEP_BETWEEN_OUTPUT_IMAGES     = 100'
    write(*,*)
  endif

  call read_value_double_precision_p(cutsnaps, 'cutsnaps')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'cutsnaps                        = 1.'
    write(*,*)
  endif

  ! jpeg images
  call read_value_logical_p(output_color_image, 'output_color_image')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'output_color_image              = .true.'
    write(*,*)
  endif

  call read_value_integer_p(imagetype_JPEG, 'imagetype_JPEG')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'imagetype_JPEG                  = 2'
    write(*,*)
  endif

  call read_value_double_precision_p(factor_subsample_image, 'factor_subsample_image')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'factor_subsample_image          = 1.0d0'
    write(*,*)
  endif

  call read_value_logical_p(USE_CONSTANT_MAX_AMPLITUDE, 'USE_CONSTANT_MAX_AMPLITUDE')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'USE_CONSTANT_MAX_AMPLITUDE      = .false.'
    write(*,*)
  endif

  call read_value_double_precision_p(CONSTANT_MAX_AMPLITUDE_TO_USE, 'CONSTANT_MAX_AMPLITUDE_TO_USE')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'CONSTANT_MAX_AMPLITUDE_TO_USE   = 1.17d4'
    write(*,*)
  endif

  call read_value_double_precision_p(POWER_DISPLAY_COLOR, 'POWER_DISPLAY_COLOR')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'POWER_DISPLAY_COLOR             = 0.30d0'
    write(*,*)
  endif

  call read_value_logical_p(DRAW_SOURCES_AND_RECEIVERS, 'DRAW_SOURCES_AND_RECEIVERS')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'DRAW_SOURCES_AND_RECEIVERS      = .true.'
    write(*,*)
  endif

  call read_value_logical_p(DRAW_WATER_IN_BLUE, 'DRAW_WATER_IN_BLUE')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'DRAW_WATER_IN_BLUE              = .true.'
    write(*,*)
  endif

  call read_value_logical_p(USE_SNAPSHOT_NUMBER_IN_FILENAME, 'USE_SNAPSHOT_NUMBER_IN_FILENAME')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'USE_SNAPSHOT_NUMBER_IN_FILENAME = .false.'
    write(*,*)
  endif

  ! postscript files
  call read_value_logical_p(output_postscript_snapshot, 'output_postscript_snapshot')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'output_postscript_snapshot      = .false.'
    write(*,*)
  endif

  call read_value_integer_p(imagetype_postscript, 'imagetype_postscript')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'imagetype_postscript            = 1'
    write(*,*)
  endif

  call read_value_logical_p(meshvect, 'meshvect')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'meshvect                        = .true.'
    write(*,*)
  endif

  call read_value_logical_p(modelvect, 'modelvect')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'modelvect                       = .false.'
    write(*,*)
  endif

  call read_value_logical_p(boundvect, 'boundvect')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'boundvect                       = .true.'
    write(*,*)
  endif

  call read_value_logical_p(interpol, 'interpol')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'interpol                        = .true.'
    write(*,*)
  endif

  call read_value_integer_p(pointsdisp, 'pointsdisp')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'pointsdisp                      = 6'
    write(*,*)
  endif

  call read_value_integer_p(subsamp_postscript, 'subsamp_postscript')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'subsamp_postscript              = 1'
    write(*,*)
  endif

  call read_value_double_precision_p(sizemax_arrows, 'sizemax_arrows')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'sizemax_arrows                  = 1.d0'
    write(*,*)
  endif

  call read_value_logical_p(US_LETTER, 'US_LETTER')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'US_LETTER                       = .false.'
    write(*,*)
  endif

  ! wavefield dumps
  call read_value_logical_p(output_wavefield_dumps, 'output_wavefield_dumps')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'output_wavefield_dumps          = .false.'
    write(*,*)
  endif

  call read_value_integer_p(imagetype_wavefield_dumps, 'imagetype_wavefield_dumps')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'imagetype_wavefield_dumps       = 1'
    write(*,*)
  endif

  call read_value_logical_p(use_binary_for_wavefield_dumps, 'use_binary_for_wavefield_dumps')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'use_binary_for_wavefield_dumps  = .false.'
    write(*,*)
  endif

  call read_value_integer_p(NUMBER_OF_SIMULTANEOUS_RUNS, 'NUMBER_OF_SIMULTANEOUS_RUNS')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'NUMBER_OF_SIMULTANEOUS_RUNS     = 1'
    write(*,*)
  endif

  call read_value_logical_p(BROADCAST_SAME_MESH_AND_MODEL, 'BROADCAST_SAME_MESH_AND_MODEL')
  if (err_occurred() /= 0) then
    some_parameters_missing_from_Par_file = .true.
    write(*,'(a)') 'BROADCAST_SAME_MESH_AND_MODEL   = .true.'
    write(*,*)
  endif

  if (some_parameters_missing_from_Par_file) then
    write(*,*)
    write(*,*) 'All the above parameters are missing from your Par_file.'
    write(*,*) 'Please cut and paste them somewhere in your Par_file (any place is fine), change their values if needed'
    write(*,*) '(the above values are just default values), and restart your run.'
    write(*,*)
    call stop_the_code('Error: some parameters are missing in your Par_file, it is incomplete or in an older format, &
       &see at the end of the standard output file of the run for detailed and easy instructions about how to fix that')
  endif

  !--------------------------------------------------------------------

  ! user output
  write(IMAIN,*) 'Title of the simulation: ',trim(title)
  write(IMAIN,*)
  if (AXISYM) write(IMAIN,*) 'Axisymmetric simulation'
  write(IMAIN,*)

  ! converts all string characters to lowercase
  irange = iachar('a') - iachar('A')
  do i = 1,len_trim(MODEL)
    if (lge(MODEL(i:i),'A') .and. lle(MODEL(i:i),'Z')) then
      MODEL(i:i) = achar(iachar(MODEL(i:i)) + irange)
    endif
  enddo
  do i = 1,len_trim(SAVE_MODEL)
    if (lge(SAVE_MODEL(i:i),'A') .and. lle(SAVE_MODEL(i:i),'Z')) then
      SAVE_MODEL(i:i) = achar(iachar(SAVE_MODEL(i:i)) + irange)
    endif
  enddo

  ! checks input parameters
  call check_parameters()

  end subroutine read_parameter_file_only

!
!-------------------------------------------------------------------------------------------------
!

  subroutine check_parameters()

  use shared_parameters

  implicit none

  ! checks partitioning
  if (NPROC <= 0) then
     print *, 'Error: Number of processes (NPROC) must be greater than or equal to one.'
     call stop_the_code('Error invalid NPROC value')
  endif

#ifndef USE_MPI
  if (NPROC > 1) then
     print *, 'Error: Number of processes (NPROC) must be equal to one when not using MPI.'
     print *, 'Please recompile with -DUSE_MPI in order to enable use of MPI.'
     call stop_the_code('Error invalid NPROC value')
  endif
#endif

  if (partitioning_method /= 1 .and. partitioning_method /= 3) then
     print *, 'Error: Invalid partitioning method number.'
     print *, 'Partitioning method ',partitioning_method,' was requested, but is not available.'
     print *, 'Support for the METIS graph partitioner has been discontinued, please use SCOTCH (option 3) instead.'
     call stop_the_code('Error invalid partitioning method')
  endif

  ! simulation parameters
  if (SIMULATION_TYPE /= 1 .and. SIMULATION_TYPE /= 3) &
    call stop_the_code('SIMULATION_TYPE can only be set to 1 or 3 in the Par_file; exiting')

  if (NOISE_TOMOGRAPHY < 0 .or. NOISE_TOMOGRAPHY > 3) &
    call stop_the_code('NOISE_TOMOGRAPHY can only be set to 0, 1, 2 or 3 in the Par_file; exiting')

  if (N_SLS < 2) call stop_the_code('must have N_SLS >= 2 even if attenuation if off because it is used to assign some arrays')

  if (ngnod /= 4 .and. ngnod /= 9) call stop_the_code('ngnod should be either 4 or 9!')

  if (seismotype < 1 .or. seismotype > 6) &
    call stop_the_code( &
'seismotype should be 1(=displ), 2(=veloc), 3(=accel), 4(=pressure), 5(=curl of displ) or 6(=the fluid potential)')

  if (USE_TRICK_FOR_BETTER_PRESSURE .and. seismotype /= 4) &
    call stop_the_code('USE_TRICK_FOR_BETTER_PRESSURE : seismograms must record pressure')

  if (subsamp_seismos < 1) call stop_the_code('Error: subsamp_seismos must be >= 1')

  if (output_color_image .and. USE_CONSTANT_MAX_AMPLITUDE .and. CONSTANT_MAX_AMPLITUDE_TO_USE < 0.d0) &
    call stop_the_code('CONSTANT_MAX_AMPLITUDE_TO_USE must be strictly positive')

  if (force_normal_to_surface .or. rec_normal_to_surface) then
    if (.not. read_external_mesh) &
      call stop_the_code('Error read_external_mesh must be set to .true. for force_normal_to_surface or rec_normal_to_surface &
            &to use external tangential_detection_curve_file')
    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) &
      call stop_the_code('NUMBER_OF_SIMULTANEOUS_RUNS not compatible with force_normal_to_surface or rec_normal_to_surface &
            &for now (look for FN2SNSR in the source code)')
  endif

  if (DT == 0.d0) call stop_the_code('DT must be non-zero value')

  ! reads in material definitions
  if (nbmodels <= 0) call stop_the_code('Non-positive number of materials not allowed!')

  ! CPML and Stacey are mutually exclusive
  if (STACEY_ABSORBING_CONDITIONS .and. PML_BOUNDARY_CONDITIONS) &
    call stop_the_code('STACEY_ABSORBING_CONDITIONS and PML_BOUNDARY_CONDITIONS are mutually exclusive but are both set to .true.')

  ! checks model
  select case (trim(MODEL))
  case ('default','ascii','binary','binary_voigt','external','gll','legacy','marmousi')
    continue ! do nothing
  case default
    print *,'Error: unknown model choosen ',trim(MODEL)
    call stop_the_code('Error bad model value for parameter MODEL')
  end select

  ! checks model
  select case (trim(SAVE_MODEL))
  case ('default','ascii','binary','gll','legacy','.false.') ! 'external' not yet
    continue ! do nothing
  case default
    print *,'Error: unknown save_model choosen ',trim(SAVE_MODEL)
    call stop_the_code('Error bad value for parameter SAVE_MODEL')
  end select

  ! check regions
  if (read_external_mesh .eqv. .false.) then
    if (nbregions <= 0) call stop_the_code('Negative number of regions not allowed for internal meshing!')
  endif

  if (NUMBER_OF_SIMULTANEOUS_RUNS <= 0) call stop_the_code('NUMBER_OF_SIMULTANEOUS_RUNS <= 0 makes no sense')

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. NPROC == 1) call stop_the_code('Serial runs require NUMBER_OF_SIMULTANEOUS_RUNS == 1')

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. trim(SAVE_MODEL) /= 'default' .and. trim(SAVE_MODEL) /= '.false.') &
    call stop_the_code('NUMBER_OF_SIMULTANEOUS_RUNS not compatible yet with SAVE_MODEL. Look for SMNSR in the source code')

  end subroutine check_parameters

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_parameter_file_receiversets()

  use constants, only: IMAIN,IIN,IN_DATA_FILES,mygroup

  use shared_parameters

  implicit none

  ! local parameters
  integer :: ireceiverlines,ier,nrec
  logical :: reread_rec_normal_to_surface
  character(len=MAX_STRING_LEN) :: stations_filename,path_to_add,dummystring

  integer,external :: err_occurred

  ! user output
  write(IMAIN,*) 'Receiver lines:'
  write(IMAIN,*) '  Nb of line sets = ',nreceiversets
  write(IMAIN,*)

  ! re-reads rec_normal_to_surface parameter to reposition read header for following next-line reads
  call read_value_logical_p(reread_rec_normal_to_surface, 'rec_normal_to_surface')
  if (err_occurred() /= 0) call stop_the_code('error reading parameter rec_normal_to_surface in Par_file')

  ! checks
  if (reread_rec_normal_to_surface .neqv. rec_normal_to_surface) call stop_the_code( &
'Invalid re-reading of rec_normal_to_surface parameter')

  ! only valid if at least 1 receiver line is specified
  if (nreceiversets < 1) call stop_the_code('number of receiver sets must be greater than 1')

  ! allocate receiver line arrays
  allocate(nrec_line(nreceiversets))
  allocate(xdeb(nreceiversets))
  allocate(zdeb(nreceiversets))
  allocate(xfin(nreceiversets))
  allocate(zfin(nreceiversets))
  allocate(record_at_surface_same_vertical(nreceiversets),stat=ier)
  if (ier /= 0 ) call stop_the_code('Error allocating receiver lines')

  nrec_line(:) = 0
  xdeb(:) = 0.d0
  zdeb(:) = 0.d0
  xfin(:) = 0.d0
  zfin(:) = 0.d0
  record_at_surface_same_vertical(:) = .false.

  ! reads in receiver sets
  if (use_existing_STATIONS) then
    write(IMAIN,*) '  using existing STATIONS file '

    ! checks if STATIONS file exisits
    stations_filename = trim(IN_DATA_FILES)//'STATIONS'

    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
      write(path_to_add,"('run',i4.4,'/')") mygroup + 1
      stations_filename = path_to_add(1:len_trim(path_to_add))//stations_filename(1:len_trim(stations_filename))
    endif

    ! counts entries
    open(unit=IIN,file=trim(stations_filename),status='old',action='read',iostat=ier)
    if (ier /= 0 ) then
      print *, 'Error could not open existing STATIONS file:'
      print *, trim(stations_filename)
      print *, 'Please check if file exists.'
      call stop_the_code('Error opening STATIONS file')
    endif
    nrec = 0
    do while(ier == 0)
      read(IIN,"(a)",iostat=ier) dummystring
      if (ier == 0) nrec = nrec + 1
    enddo
    close(IIN)

    write(IMAIN,*) '  file name is ',trim(stations_filename)
    write(IMAIN,*) '  found ',nrec,' receivers'
    write(IMAIN,*)

  else
    ! loop on all the receiver lines
    do ireceiverlines = 1,nreceiversets
      call read_value_integer_next_p(nrec_line(ireceiverlines),'nrec')
      if (err_occurred() /= 0) call stop_the_code('error reading parameter nrec in Par_file')

      call read_value_double_prec_next_p(xdeb(ireceiverlines),'xdeb')
      if (err_occurred() /= 0) call stop_the_code('error reading parameter xdeb in Par_file')

      call read_value_double_prec_next_p(zdeb(ireceiverlines),'zdeb')
      if (err_occurred() /= 0) call stop_the_code('error reading parameter zdeb in Par_file')

      call read_value_double_prec_next_p(xfin(ireceiverlines),'xfin')
      if (err_occurred() /= 0) call stop_the_code('error reading parameter xfin in Par_file')

      call read_value_double_prec_next_p(zfin(ireceiverlines),'zfin')
      if (err_occurred() /= 0) call stop_the_code('error reading parameter zfin in Par_file')

      call read_value_logical_next_p(record_at_surface_same_vertical(ireceiverlines),'record_at_surface_same_vertical')
      if (err_occurred() /= 0) call stop_the_code('error reading parameter record_at_surface_same_vertical in Par_file')

      if (read_external_mesh .and. record_at_surface_same_vertical(ireceiverlines)) then
        call stop_the_code('Cannot use record_at_surface_same_vertical with external meshes!')
      endif
    enddo
  endif

  end subroutine read_parameter_file_receiversets

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_parameter_file_derive_flags()

  use shared_parameters

  implicit none

  ! derives additional flags based on input parameters

  ! sets overal Bielak flag
  add_Bielak_conditions = add_Bielak_conditions_bottom .or. add_Bielak_conditions_right .or. &
                          add_Bielak_conditions_top .or. add_Bielak_conditions_left

  ! boundary conditions
  if (add_Bielak_conditions .and. .not. STACEY_ABSORBING_CONDITIONS) &
    call stop_the_code('need STACEY_ABSORBING_CONDITIONS set to .true. in order to use add_Bielak_conditions')

  ! solve the conflict in value of PML_BOUNDARY_CONDITIONS and STACEY_ABSORBING_CONDITIONS
  if (PML_BOUNDARY_CONDITIONS) any_abs = .true.
  if (STACEY_ABSORBING_CONDITIONS) any_abs = .true.

  ! initializes flags for absorbing boundaries
  if (.not. any_abs) then
    absorbbottom = .false.
    absorbright = .false.
    absorbtop = .false.
    absorbleft = .false.
  endif

  ! can use only one point to display lower-left corner only for interpolated snapshot
  if (pointsdisp < 3) then
    pointsdisp = 3
    plot_lowerleft_corner_only = .true.
  else
    plot_lowerleft_corner_only = .false.
  endif

  end subroutine read_parameter_file_derive_flags
