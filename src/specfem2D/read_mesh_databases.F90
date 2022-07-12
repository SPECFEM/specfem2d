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

  subroutine read_mesh_for_init()

  use constants, only: &
    IMAIN,IIN,DISPLAY_COLORS,DISPLAY_ELEMENT_NUMBERS_POSTSCRIPT,OUTPUT_FILES, &
    ADD_RANDOM_PERTURBATION_TO_THE_MESH,ADD_PERTURBATION_AROUND_SOURCE_ONLY

  use specfem_par
  use specfem_par_movie

  implicit none

  ! local parameters
  integer, external :: err_occurred
  integer :: npgeo_all
  integer :: ier,int_dummy
  character(len=MAX_STRING_LEN) :: prname, dummy

  logical :: do_rerun_mesher
  integer :: local_i
  integer (kind=RegInt_K) :: local_ireg
  double precision :: local_dble
  logical :: local_l,local_l1,local_l2,local_l3
  character(len=MAX_STRING_LEN) :: local_str

  ! starts reading SIMULATION_TYPE and SAVE_FORWARD from Par_file

  call open_parameter_file()

  call read_value_string_p(dummy, 'solver.title')
  if (err_occurred() /= 0) call stop_the_code('error reading parameter title in Par_file')

  ! read type of simulation
  call read_value_integer_p(SIMULATION_TYPE, 'solver.SIMULATION_TYPE')
  if (err_occurred() /= 0) call stop_the_code('error reading parameter SIMULATION_TYPE in Par_file')

  call read_value_integer_p(int_dummy, 'solver.NOISE_TOMOGRAPHY')
  if (err_occurred() /= 0) call stop_the_code('error reading parameter NOISE_TOMOGRAPHY in Par_file')

  call read_value_logical_p(SAVE_FORWARD, 'solver.SAVE_FORWARD')
  if (err_occurred() /= 0) call stop_the_code('error reading parameter SAVE_FORWARD in Par_file')

  call close_parameter_file()

  ! starts reading in parameters from input Database file

  ! opens Database file
  write(prname,"(a,i5.5,a)") trim(OUTPUT_FILES)//'Database',myrank,'.bin'
  ! note: adding access='stream' would further decrease file size
  open(unit=IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0 ) then
    if (myrank == 0) then
      print *,'Error opening database file: ',trim(prname)
      print *
      print *,'Please make sure that the mesher has been run before this solver simulation with the correct settings...'
    endif
    call exit_MPI(myrank,'Error opening file '//trim(OUTPUT_FILES)//'Database***.bin')
  endif

! note: we now read in the parameter file also in the initialization of the solver.
!       the idea is to separate solver and mesher further and
!       to be able to run only the solver if Par_file parameters were changed which don't affect the meshing stage.
!       this should include parameters like simulation_type, sources, receivers, image outputs etc. and would be more
!       similar of how to run simulations in the 3D version.
!
!       however, at the moment most of the Par_file changes are setup in the mesher and require re-running the mesher.
!
!       we will use local_** variables just to read in the settings stored in the database and then compare with
!       the ones from the Par_file to decide if we need to re-run the mesher.
  do_rerun_mesher = .false.

  simulation_title = trim(title) ! from parameter file

  !-------- starts reading init section

  !---  read job title and skip remaining titles of the input file
  read(IIN) local_str ! simulation_title
  !print *,'debug: simulation_title',myrank,'***',local_str,'***'

  !---- read parameters from input file
  read(IIN) local_i, local_l ! NOISE_TOMOGRAPHY, UNDO_ATTENUATION_AND_OR_PML
  !print *,'debug: NOISE_TOMOGRAPHY,UNDO_ATTENUATION_AND_OR_PML ',myrank,local_i,local_l

  read(IIN) nspec
  !print *,'debug: nspec ',myrank,nspec

  read(IIN) npgeo,nproc_read_from_database
  !print *,'debug: npgeo,nproc_read_from_database ',myrank,npgeo,nproc_read_from_database

  read(IIN) local_l1, local_l2 ! output_grid_Gnuplot,interpol
  !print *,'debug: output_grid_Gnuplot,interpol',myrank,local_l1,local_l2

  read(IIN) local_i ! NTSTEP_BETWEEN_OUTPUT_INFO
  !print *,'debug: NTSTEP_BETWEEN_OUTPUT_INFO ',myrank,local_i

  read(IIN) local_i ! NTSTEP_BETWEEN_OUTPUT_SEISMOS
  !print *,'debug: NTSTEP_BETWEEN_OUTPUT_SEISMOS ',myrank,local_i

  read(IIN) local_i ! NTSTEP_BETWEEN_OUTPUT_IMAGES
  !print *,'debug: NTSTEP_BETWEEN_OUTPUT_IMAGES ',myrank,local_i

  read(IIN) local_l ! PML_BOUNDARY_CONDITIONS
  !print *,'debug: PML_BOUNDARY_CONDITIONS ',myrank,local_l

  if (local_l .neqv. PML_BOUNDARY_CONDITIONS) then
    print *,'Warning: rank ',myrank,' read mesh: PML_BOUNDARY_CONDITIONS setting changed'
    do_rerun_mesher = .true.
  endif

  read(IIN) local_l ! ROTATE_PML_ACTIVATE

  read(IIN) local_dble ! ROTATE_PML_ANGLE

  read(IIN) local_dble ! K_MIN_PML

  read(IIN) local_dble ! K_MAX_PML

  read(IIN) local_dble ! damping_change_factor_acoustic

  read(IIN) local_dble ! damping_change_factor_elastic

  read(IIN) local_l ! PML_PARAMETER_ADJUSTMENT

  read(IIN) local_l ! read_external_mesh
  if (local_l .neqv. read_external_mesh) then
     print *,'Warning: rank ',myrank,' read mesh: read_external_mesh setting changed'
     do_rerun_mesher = .true.
   endif

  read(IIN) local_i ! NELEM_PML_THICKNESS
  if (local_i /= NELEM_PML_THICKNESS) then
     print *,'Warning: rank ',myrank,' read mesh: NELEM_PML_THICKNESS setting changed'
     do_rerun_mesher = .true.
   endif

  read(IIN) local_i, local_i, local_i ! NTSTEP_BETWEEN_OUTPUT_SAMPLE,imagetype_JPEG,imagetype_wavefield_dumps

  read(IIN) local_l, local_l ! output_postscript_snapshot,output_color_image

  ! meshvect,modelvect,boundvect,cutsnaps,subsamp_postscript,sizemax_arrows
  read(IIN) local_l, local_l, local_l, local_dble, local_i, local_dble

  read(IIN) local_dble ! anglerec

  read(IIN) local_l ! initialfield
  if (local_l .neqv. initialfield) then
     print *,'Warning: rank ',myrank,' read mesh: initialfield setting changed'
     do_rerun_mesher = .true.
   endif

  ! add_Bielak_conditions_bottom,add_Bielak_conditions_right,add_Bielak_conditions_top,add_Bielak_conditions_left
  read(IIN) local_l,local_l,local_l,local_l

  read(IIN) local_str,local_i ! seismotype,imagetype_postscript

  read(IIN) local_str ! MODEL
  if (trim(local_str) /= trim(MODEL)) then
     print *,'Warning: rank ',myrank,' read mesh: MODEL setting changed'
     do_rerun_mesher = .true.
   endif

  read(IIN) local_str ! SAVE_MODEL

  read(IIN) local_str ! TOMOGRAPHY_FILE
  if (trim(local_str) /= trim(TOMOGRAPHY_FILE)) then
     print *,'Warning: rank ',myrank,' read mesh: TOMOGRAPHY_FILE setting changed'
     do_rerun_mesher = .true.
   endif

  ! output_grid_ASCII,OUTPUT_ENERGY,NTSTEP_BETWEEN_OUTPUT_ENERGY,output_wavefield_dumps
  read(IIN) local_l, local_l, local_i, local_l

  read(IIN) local_l ! use_binary_for_wavefield_dumps

  read(IIN) local_l1, local_l2, local_l3 ! ATTENUATION_VISCOELASTIC,ATTENUATION_PORO_FLUID_PART,ATTENUATION_VISCOACOUSTIC
  if (local_l1 .neqv. ATTENUATION_VISCOELASTIC) then
     print *,'Warning: rank ',myrank,' read mesh: ATTENUATION_VISCOELASTIC setting changed'
     do_rerun_mesher = .true.
   endif       ! elements load changes in meshing
  if (local_l2 .neqv. ATTENUATION_PORO_FLUID_PART) then
     print *,'Warning: rank ',myrank,' read mesh: ATTENUATION_PORO_FLUID_PART setting changed'
     do_rerun_mesher = .true.
   endif
  if (local_l3 .neqv. ATTENUATION_VISCOACOUSTIC) then
     print *,'Warning: rank ',myrank,' read mesh: ATTENUATION_VISCOACOUSTIC setting changed'
     do_rerun_mesher = .true.
   endif

  read(IIN) local_l ! USE_SOLVOPT

  read(IIN) local_l ! save_ASCII_seismograms

  read(IIN) local_l, local_l ! save_binary_seismograms_single,save_binary_seismograms_double

  read(IIN) local_l ! USE_TRICK_FOR_BETTER_PRESSURE

  read(IIN) local_l ! COMPUTE_INTEGRATED_ENERGY_FIELD

  read(IIN) local_l ! save_ASCII_kernels

  read(IIN) local_i ! NTSTEP_BETWEEN_COMPUTE_KERNELS

  read(IIN) local_l ! APPROXIMATE_HESS_KL

  read(IIN) local_l ! NO_BACKWARD_RECONSTRUCTION

  read(IIN) local_l ! DRAW_SOURCES_AND_RECEIVERS

  read(IIN) local_dble,local_dble ! Q0_poroelastic,freq0_poroelastic

  read(IIN) local_l ! AXISYM
  if (local_l .neqv. AXISYM) then
     print *,'Warning: rank ',myrank,' read mesh: AXISYM setting changed'
     do_rerun_mesher = .true.
   endif

  read(IIN) local_l ! P_SV

  read(IIN) local_dble ! factor_subsample_image

  read(IIN) local_l ! USE_CONSTANT_MAX_AMPLITUDE

  read(IIN) local_dble ! CONSTANT_MAX_AMPLITUDE_TO_USE

  read(IIN) local_l ! USE_SNAPSHOT_NUMBER_IN_FILENAME

  read(IIN) local_l ! DRAW_WATER_IN_BLUE

  read(IIN) local_l ! US_LETTER

  read(IIN) local_dble ! POWER_DISPLAY_COLOR

  read(IIN) local_l ! SU_FORMAT

  read(IIN) local_dble ! USER_T0

  read(IIN) local_i ! time_stepping_scheme

  read(IIN) local_l ! ADD_PERIODIC_CONDITIONS
  if (local_l .neqv. ADD_PERIODIC_CONDITIONS) then
     print *,'Warning: rank ',myrank,' read mesh: ADD_PERIODIC_CONDITIONS setting changed'
     do_rerun_mesher = .true.
   endif

  read(IIN) local_dble ! PERIODIC_HORIZ_DIST
  if (local_dble /= PERIODIC_HORIZ_DIST) then
     print *,'Warning: rank ',myrank,' read mesh: PERIODIC_HORIZ_DIST setting changed'
     do_rerun_mesher = .true.
   endif

  read(IIN) local_l ! GPU_MODE

  read(IIN) local_i ! setup_with_binary_database

  !---- read time step
  read(IIN) local_ireg,local_dble ! NSTEP,DT

  read(IIN) local_i ! NT_DUMP_ATTENUATION

  ! read the ACOUSTIC_FORCING flag
  read(IIN) local_l ! ACOUSTIC_FORCING
  if (local_l .neqv. ACOUSTIC_FORCING) then
     print *,'Warning: rank ',myrank,' read mesh: ACOUSTIC_FORCING setting changed'
     do_rerun_mesher = .true.
   endif

  ! 'NUMBER_OF_SIMULTANEOUS_RUNS'
  read(IIN) local_i ! NUMBER_OF_SIMULTANEOUS_RUNS

  ! 'BROADCAST_SAME_MESH_AND_MODEL'
  read(IIN) local_l ! BROADCAST_SAME_MESH_AND_MODEL

  read(IIN) local_l1,local_l2 ! ADD_RANDOM_PERTURBATION_TO_THE_MESH,ADD_PERTURBATION_AROUND_SOURCE_ONLY
  if ((local_l1 .neqv. ADD_RANDOM_PERTURBATION_TO_THE_MESH) &
      .or. (local_l2 .neqv. ADD_PERTURBATION_AROUND_SOURCE_ONLY)) then
    print *,'Warning: rank ',myrank,' read mesh: ADD_RANDOM_PERTURBATION or ADD_PERTURBATION_AROUND_SOURCE setting changed'
    do_rerun_mesher = .true.
  endif

  ! check if setup needs re-running mesher
  if (do_rerun_mesher) then
    write(IMAIN,*)
    write(IMAIN,*) 'Setup in Par_file changed and requires re-running the mesher...exiting'
    write(IMAIN,*)
    call stop_the_code('Please re-run mesher')
  endif
  call synchronize_all()

  ! At that point seismotype is a string (e.g 2,3,6). The following routine convert it to an integer array: seismotypeVec
  call process_seismotype_line()

  !-------- finish reading init section
  ! collect infos from all
  call sum_all_i(npgeo,npgeo_all)

  ! sets time step for time scheme
  deltat = real(DT,kind=CUSTOM_REAL)

  ! user output
  if (myrank == 0) then
    !---- print the date, time and start-up banner
    call datim(simulation_title)

    ! AXISYM simulation
    if (AXISYM) then
      write(IMAIN,*)
      write(IMAIN,*)
      write(IMAIN,*) '-----------------------------------------------------'
      write(IMAIN,*) '--- A x i s y m m e t r i c   S i m u l a t i o n ---'
      write(IMAIN,*) '-----------------------------------------------------'
      call flush_IMAIN()
    endif

    ! outputs parameters read
    write(IMAIN,200) npgeo_all,NDIM
    write(IMAIN,600) NTSTEP_BETWEEN_OUTPUT_INFO,DISPLAY_COLORS,DISPLAY_ELEMENT_NUMBERS_POSTSCRIPT
    write(IMAIN,700) trim(seismotype),anglerec
    write(IMAIN,750) initialfield, &
                 add_Bielak_conditions_bottom,add_Bielak_conditions_right,add_Bielak_conditions_top,add_Bielak_conditions_left, &
                 ATTENUATION_VISCOELASTIC,ATTENUATION_VISCOACOUSTIC,output_grid_ASCII,OUTPUT_ENERGY
    write(IMAIN,800) imagetype_postscript,cutsnaps,subsamp_postscript

    ! time step
    write(IMAIN,703) NSTEP,deltat,NSTEP*deltat

    ! flush buffer
    call flush_IMAIN()
  endif

  ! -------------- formatting
  ! output formats
200 format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Number of spectral element control nodes. . .(npgeo) =',i8/5x, &
  'Number of space dimensions. . . . . . . . . . (NDIM) =',i8)

600 format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Display frequency . . . (NTSTEP_BETWEEN_OUTPUT_INFO) = ',i6/ 5x, &
  'Color display . . . . . . . . . . . . . . . (colors) = ',i6/ 5x, &
  ' == 0     black and white display              ',  / 5x, &
  ' == 1     color display                        ',  /5x, &
  'Numbered mesh . . . . . . . . . . . . . . .(numbers) = ',i6/ 5x, &
  ' == 0     do not number the mesh               ',  /5x, &
  ' == 1     number the mesh                      ')

700 format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Seismograms recording type . . . . . . .(seismotype) = ',a/5x, &
  'Angle for first line of receivers. . . . .(anglerec) = ',f6.2)

750 format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Read external initial field. . . . . .(initialfield) = ',l6/5x, &
  'Add Bielak conditions (add_Bielak_conditions_bottom) = ',l6/5x, &
  'Add Bielak conditions .(add_Bielak_conditions_right) = ',l6/5x, &
  'Add Bielak conditions . .(add_Bielak_conditions_top) = ',l6/5x, &
  'Add Bielak conditions. .(add_Bielak_conditions_left) = ',l6/5x, &
  'Attenuation in solid . . .(ATTENUATION_VISCOELASTIC) = ',l6/5x, &
  'Attenuation in fluid . . (ATTENUATION_VISCOACOUSTIC) = ',l6/5x, &
  'Save grid in ASCII file or not . (output_grid_ASCII) = ',l6/5x, &
  'Save a file with total energy or not.(OUTPUT_ENERGY) = ',l6)

800 format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Vector display type . . . . . . . . .(imagetype_postscript) = ',i6/5x, &
  'Percentage of cut for vector plots. . . . . . . .(cutsnaps) = ',f6.2/5x, &
  'Subsampling of velocity model display. (subsamp_postscript) = ',i6)

703 format(//' I t e r a t i o n s '/1x,19('='),//5x, &
      'Number of time iterations . . . . .(NSTEP) =',i8,/5x, &
      'Time step increment. . . . . . . . . .(DT) =',1pe15.6,/5x, &
      'Total simulation duration . . . . . (ttot) =',1pe15.6)

  end subroutine read_mesh_for_init


!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases()

  use constants, only: IIN,IMAIN,ADD_A_SMALL_CRACK_IN_THE_MEDIUM
  use specfem_par

  implicit none

  ! continues reading database files
  ! note: we opened the database file in read_mesh_for_init() and will continue reading from its current position
  !       thus, the ordering here must be consistent with the order in save_databases.f90

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*)
    write(IMAIN,*) 'reading mesh databases:'
    call flush_IMAIN()
  endif

  ! reads the spectral macrobloc nodal coordinates
  ! and basic properties of the spectral elements
  call read_mesh_databases_coorg_elem()

  ! reads attenuation information (needs source timefunction type from read sources before)
  call read_mesh_databases_attenuation()

  ! material properties
  call read_mesh_databases_mato()

  ! add a small crack (discontinuity) in the medium manually
  if (ADD_A_SMALL_CRACK_IN_THE_MEDIUM) call add_manual_crack()

  ! determines if each spectral element is elastic, poroelastic, or acoustic
  call get_simulation_domains()

  ! reads interfaces data
  call read_mesh_databases_interfaces()

  ! reads absorbing boundary data
  call read_mesh_databases_absorbing()

  ! reads acoustic forcing boundary data
  call read_mesh_databases_acoustic_forcing()

  ! reads acoustic free surface data
  call read_mesh_databases_free_surf()

  ! reads coupled edges
  call read_mesh_databases_coupled()

  ! reads tangential detection curve
  call read_mesh_databases_tangential()

  ! reads axial elements data
  call read_mesh_databases_axial_elements()

  ! end of reading
  ! closes input Database file
  close(IIN)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'done reading mesh databases'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine read_mesh_databases


!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_coorg_elem()

! reads the spectral macrobloc nodal coordinates

  use constants, only: IMAIN,IIN
  use specfem_par
  use specfem_par_movie

  implicit none

  ! local parameters
  integer :: ipoin,ip,ier
  double precision, dimension(NDIM) :: coorgread
  integer :: nspec_all,nelem_acforcing_all,nelem_acoustic_surface_all

  ! user output
  if (myrank == 0) write(IMAIN,*) 'reading nodal coordinates...'

  ! safety check
  if (NDIM /= 2) call stop_the_code('Invalid NDIM value to read coordinates, please check...')

  ! allocates nodal coordinates
  allocate(coorg(NDIM,npgeo),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating coorg array')

  ! initializes
  coorg(:,:) = 0.d0

  ! reads the spectral macrobloc nodal coordinates

  ! reads in values
  ipoin = 0
  do ip = 1,npgeo
    ! reads coordinates
    read(IIN) ipoin,coorgread(1),coorgread(2)

    ! checks index
    if (ipoin < 1 .or. ipoin > npgeo) then
      print *, 'Error reading coordinates: invalid point number',ipoin,coorgread(1),coorgread(2), &
               'at position',ip,'out of',npgeo
      call exit_MPI(myrank,'Wrong control point number')
    endif

    ! saves coordinate array
    coorg(:,ipoin) = coorgread(:)
  enddo

  !---- read the basic properties of the spectral elements
  read(IIN) numat,NGNOD,nspec,pointsdisp,plot_lowerleft_corner_only

  read(IIN) nelemabs,nelem_acforcing,nelem_acoustic_surface,num_fluid_solid_edges, &
              num_fluid_poro_edges,num_solid_poro_edges,nnodes_tangential_curve, &
              nelem_on_the_axis

  ! collects numbers
  call sum_all_i(nspec,nspec_all)
  call sum_all_i(nelem_acforcing,nelem_acforcing_all)
  call sum_all_i(nelem_acoustic_surface,nelem_acoustic_surface_all)

  ! user output
  if (myrank == 0) then
    ! print element group main parameters
    write(IMAIN,107)
    write(IMAIN,207) nspec_all,NGNOD,NGLLX,NGLLZ,NGLLX*NGLLZ,pointsdisp,numat, &
                     nelem_acforcing_all,nelem_acoustic_surface_all
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! allocates mesh arrays
  ! default material
  allocate(density(2,numat))
  allocate(anisotropycoef(10,numat)) ! don't forget c22 value (it is used for AXISYM simulations only)

  allocate(porosity(numat), &
           tortuosity(numat), &
           permeability(3,numat),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating porosity,.. arrays')

  allocate(poroelastcoef(4,3,numat),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating poroelastcoef arrays')

  allocate(QKappa_attenuationcoef(numat), &
           Qmu_attenuationcoef(numat),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating attenuation arrays')


  ! output formats
107 format(/1x,'-- Spectral Elements --',//)

207 format(5x,'Number of spectral elements . . . . . . . . .  (nspec) =',i7,/5x, &
               'Number of control nodes per element . . . . . (NGNOD) =',i7,/5x, &
               'Number of points in X-direction . . . . . . . (NGLLX) =',i7,/5x, &
               'Number of points in Y-direction . . . . . . . (NGLLZ) =',i7,/5x, &
               'Number of points per element. . . . . . (NGLLX*NGLLZ) =',i7,/5x, &
               'Number of points for display . . . . . . (pointsdisp) =',i7,/5x, &
               'Number of element material sets . . . . . . . (numat) =',i7,/5x, &
               'Number of acoustic forcing elements (nelem_acforcing) =',i7,/5x, &
               'Number of acoustic free surf (nelem_acoustic_surface) =',i7)

  end subroutine read_mesh_databases_coorg_elem

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_attenuation()

! reads attenuation information

  use constants, only: IIN,IMAIN,NGLLX,NGLLZ

  use specfem_par, only: myrank,N_SLS,ATTENUATION_f0_REFERENCE,READ_VELOCITIES_AT_f0

  implicit none

  ! user output
  if (myrank == 0) write(IMAIN,*) 'reading attenuation setup...'

  ! attenuation parameters
  read(IIN) N_SLS, ATTENUATION_f0_REFERENCE, READ_VELOCITIES_AT_f0

  ! checks number of standard linear solids
  if (N_SLS < 1) call stop_the_code('must have N_SLS >= 1 even if attenuation if off because it is used to assign some arrays')

  end subroutine read_mesh_databases_attenuation

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_mato()

! reads spectral macrobloc data

  use constants, only: IIN,IMAIN,myrank

  use specfem_par, only: nspec,NGNOD,kmato,knods,region_CPML,f0_source,numat

  implicit none

  ! local parameters
  integer :: n,k,ispec,kmato_read,pml_read,imat,ier
  integer, dimension(:), allocatable :: knods_read

  ! user output
  if (myrank == 0) write(IMAIN,*) 'reading material properties..'

  ! reads and sets the material properties
  call read_materials(f0_source(1))

  ! add support for using PML in MPI mode with external mesh
  allocate(region_CPML(nspec),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating region_CPML array')

  ! elements
  allocate(kmato(nspec), &
           knods(NGNOD,nspec),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating kmato,.. array')

  ! initializes
  kmato(:) = 0
  knods(:,:) = 0
  region_CPML(:) = 0

  ! temporary read array
  allocate(knods_read(NGNOD),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating temporary knods array')

  ! reads spectral macrobloc data

  ! reads in values
  n = 0
  do ispec = 1,nspec
    ! format: #element_id  #material_id #node_id1 #node_id2 #...
    read(IIN) n,kmato_read,(knods_read(k),k = 1,NGNOD),pml_read

    ! material association
    kmato(n) = kmato_read
    region_CPML(n) = pml_read

    ! element control node indices
    knods(:,n)= knods_read(:)
  enddo

  ! checks material array
  do ispec = 1,nspec
    imat = kmato(ispec)
    if (imat < 1 .or. imat > numat) then
      print *,'Error: rank ',myrank,'found element ',ispec,' with invalid material id',imat
      print *,'Please check your material definitions and element assignments'
      call stop_the_code('Invalid material id found')
    endif
  enddo

  ! frees temporary array
  deallocate(knods_read)

  end subroutine read_mesh_databases_mato

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_interfaces()

! reads in interface dimensions

  use constants, only: IIN,IMAIN
  use specfem_par

  implicit none

  ! local parameters
  integer :: num_interface,ie,ier,my_interfaces_read

  ! user output
  if (myrank == 0) write(IMAIN,*) 'reading interfaces informations...'

  ! reads number of interfaces
  read(IIN) ninterface, max_interface_size

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  number of interfaces         = ',ninterface
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! allocates arrays for mpi/partition interfaces
  if (ninterface > 0) then
    allocate(my_neighbors(ninterface), &
             my_nelmnts_neighbors(ninterface), &
             my_interfaces(4,max_interface_size,ninterface),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating interfaces arrays')

    allocate(ibool_interfaces_acoustic(NGLLX*max_interface_size,ninterface), &
             ibool_interfaces_elastic(NGLLX*max_interface_size,ninterface), &
             ibool_interfaces_poroelastic(NGLLX*max_interface_size,ninterface), &
             ibool_interfaces_ext_mesh_init(NGLLX*max_interface_size,ninterface),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating interfaces ibool arrays')

    allocate(nibool_interfaces_acoustic(ninterface), &
             nibool_interfaces_elastic(ninterface), &
             nibool_interfaces_poroelastic(ninterface), &
             nibool_interfaces_ext_mesh(ninterface), &
             inum_interfaces_acoustic(ninterface), &
             inum_interfaces_elastic(ninterface), &
             inum_interfaces_poroelastic(ninterface),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating interfaces nibool arrays')
  else
    ! dummy
    allocate(my_neighbors(1),my_nelmnts_neighbors(1),my_interfaces(1,1,1))
    allocate(ibool_interfaces_acoustic(1,1),ibool_interfaces_elastic(1,1))
    allocate(ibool_interfaces_poroelastic(1,1),ibool_interfaces_ext_mesh_init(1,1))
    allocate(nibool_interfaces_acoustic(1),nibool_interfaces_elastic(1))
    allocate(nibool_interfaces_poroelastic(1),nibool_interfaces_ext_mesh(1))
    allocate(inum_interfaces_acoustic(1),inum_interfaces_elastic(1),inum_interfaces_poroelastic(1))
  endif
  ! initializes
  my_neighbors(:) = -1
  my_nelmnts_neighbors(:) = 0
  my_interfaces(:,:,:) = -1

  ! note: for serial simulations, ninterface will be zero.
  !       thus no further reading will be done below

  ! reads in interfaces
  do num_interface = 1, ninterface
    ! format: #process_interface_id  #number_of_elements_on_interface
    ! where
    !     process_interface_id = rank of (neighbor) process to share MPI interface with
    !     number_of_elements_on_interface = number of interface elements
    read(IIN) my_neighbors(num_interface), my_nelmnts_neighbors(num_interface)

    ! loops over interface elements
    do ie = 1, my_nelmnts_neighbors(num_interface)
      ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2
      !
      ! interface types:
      !     1  -  corner point only
      !     2  -  element edge
      read(IIN) my_interfaces_read, my_interfaces(2,ie,num_interface), &
                my_interfaces(3,ie,num_interface), my_interfaces(4,ie,num_interface)

      my_interfaces(1,ie,num_interface) = my_interfaces_read
    enddo
  enddo

  end subroutine read_mesh_databases_interfaces


!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_absorbing()

! reads in absorbing edges

  use constants, only: IMAIN,IIN,IEDGE1,IEDGE2,IEDGE3,IEDGE4
  use specfem_par

  implicit none

  ! local parameters
  integer :: inum,inum_duplicate,numabsread,typeabsread
  logical :: codeabsread(4)

  integer :: num_abs_boundary_faces_tot,nspec_left_tot,nspec_right_tot,nspec_bottom_tot,nspec_top_tot
  integer :: ier,num_all
  integer, dimension(:), allocatable :: numabs

  ! user output
  if (myrank == 0) write(IMAIN,*) 'reading absorbing boundary...'

  ! saftey check
  if (nelemabs < 0) then
    if (myrank == 0) write(IMAIN,*) '  Warning: read in negative nelemabs ',nelemabs,'...resetting to zero!'
    nelemabs = 0
  endif

  ! number of boundary element faces (uses same name as in 3D versions; in 2D, the boundary would be edges)
  ! (making sure that num_abs_boundary_faces is zero if there are no boundary faces)
  num_abs_boundary_faces = nelemabs

  ! gets total number of boundary faces/edges
  call sum_all_i(num_abs_boundary_faces,num_all)

  ! determines flag for absorbing boundaries
  if (num_abs_boundary_faces <= 0) then
    anyabs = .false.
  else
    anyabs = .true.
  endif

  ! note: below is commented out as it would overwrite the setting coming from Par_file and make this setting
  !       local, i.e., different for different MPI slices. then, any MPI communication within an if-statement could
  !       potentially stall the execution. thus, it's a dangerous handling. Instead, we could add a new, local flag
  !       to make it more explicit when/when not to handle stacey routine calls - left as a possible todo in future.
  !
  ! sets (local) Stacey flag
  !if (anyabs .and. (.not. PML_BOUNDARY_CONDITIONS)) then
  !  STACEY_ABSORBING_CONDITIONS = .true.
  !else
  !  STACEY_ABSORBING_CONDITIONS = .false.
  !endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  Total number of absorbing elements = ',num_all
    write(IMAIN,*)
    write(IMAIN,*) '  any absorbing boundary flag        = ',anyabs
    write(IMAIN,*) '  PML boundary flag                  = ',PML_BOUNDARY_CONDITIONS
    write(IMAIN,*) '  Stacey boundary flag               = ',STACEY_ABSORBING_CONDITIONS
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! temporary array
  if (num_abs_boundary_faces > 0) then
    allocate(numabs(num_abs_boundary_faces),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating temporary absorbing array')
  else
    ! dummy
    allocate(numabs(1))
  endif
  numabs(:) = 0

  ! allocate arrays for absorbing boundary conditions
  if (num_abs_boundary_faces > 0) then
    allocate(codeabs(4,num_abs_boundary_faces),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating absorbing array')
  else
    ! dummy
    allocate(codeabs(1,1))
  endif
  codeabs(:,:) = .false.

  !---codeabs_corner(1,num_abs_boundary_faces) denotes whether element is on bottom-left corner of absorbing boundary or not
  !---codeabs_corner(2,num_abs_boundary_faces) denotes whether element is on bottom-right corner of absorbing boundary or not
  !---codeabs_corner(3,num_abs_boundary_faces) denotes whether element is on top-left corner of absorbing boundary or not
  !---codeabs_corner(4,num_abs_boundary_faces) denotes whether element is on top-right corner of absorbing boundary or not
  if (num_abs_boundary_faces > 0) then
    allocate(codeabs_corner(4,num_abs_boundary_faces),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating absorbing codeabs_corner array')
  else
    ! dummy
    allocate(codeabs_corner(1,1))
  endif
  codeabs_corner(:,:) = .false.

  if (num_abs_boundary_faces > 0) then
    allocate(abs_boundary_type(num_abs_boundary_faces), &
             ibegin_edge1(num_abs_boundary_faces), &
             iend_edge1(num_abs_boundary_faces), &
             ibegin_edge2(num_abs_boundary_faces), &
             iend_edge2(num_abs_boundary_faces), &
             ibegin_edge3(num_abs_boundary_faces), &
             iend_edge3(num_abs_boundary_faces), &
             ibegin_edge4(num_abs_boundary_faces), &
             iend_edge4(num_abs_boundary_faces),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating absorbing boundary index arrays')
    ! poroelastic
    allocate(ibegin_edge1_poro(num_abs_boundary_faces), &
             iend_edge1_poro(num_abs_boundary_faces), &
             ibegin_edge2_poro(num_abs_boundary_faces), &
             iend_edge2_poro(num_abs_boundary_faces), &
             ibegin_edge3_poro(num_abs_boundary_faces), &
             iend_edge3_poro(num_abs_boundary_faces), &
             ibegin_edge4_poro(num_abs_boundary_faces), &
             iend_edge4_poro(num_abs_boundary_faces),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating absorbing boundary poro arrays')

    allocate(ib_left(num_abs_boundary_faces), &
             ib_right(num_abs_boundary_faces), &
             ib_bottom(num_abs_boundary_faces), &
             ib_top(num_abs_boundary_faces),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating absorbing boundary arrays')
  else
    ! dummy
    allocate(abs_boundary_type(1),ibegin_edge1(1),iend_edge1(1),ibegin_edge2(1),iend_edge2(1))
    allocate(ibegin_edge3(1),iend_edge3(1),ibegin_edge4(1),iend_edge4(1))
    allocate(ibegin_edge1_poro(1),iend_edge1_poro(1),ibegin_edge2_poro(1),iend_edge2_poro(1))
    allocate(ibegin_edge3_poro(1),iend_edge3_poro(1),ibegin_edge4_poro(1),iend_edge4_poro(1))
    allocate(ib_left(1),ib_right(1),ib_bottom(1),ib_top(1))
  endif

  ! initializes
  abs_boundary_type(:) = 0

  ibegin_edge1(:) = 0
  iend_edge1(:) = 0

  ibegin_edge2(:) = 0
  iend_edge2(:) = 0

  ibegin_edge3(:) = 0
  iend_edge3(:) = 0

  ibegin_edge4(:) = 0
  iend_edge4(:) = 0

  ! poroelastic edges (have not been set in mesher yet)
  ibegin_edge1_poro(:) = 1
  iend_edge1_poro(:) = NGLLX

  ibegin_edge2_poro(:) = 1
  iend_edge2_poro(:) = NGLLZ

  ibegin_edge3_poro(:) = 1
  iend_edge3_poro(:) = NGLLX

  ibegin_edge4_poro(:) = 1
  iend_edge4_poro(:) = NGLLZ

  nspec_left = 0
  nspec_right = 0
  nspec_bottom = 0
  nspec_top = 0

  ib_right(:) = 0
  ib_left(:) = 0
  ib_bottom(:) = 0
  ib_top(:) = 0

  ! reads in absorbing edges
  if (anyabs) then

    ! reads absorbing boundaries
    do inum = 1,num_abs_boundary_faces

      ! beware here and below that external meshes (for instance coming from CUBIT or Gmsh)
      ! may have rotated elements and thus edge 1 may not correspond to the bottom,
      ! edge 2 may not correspond to the right, edge 3 may not correspond to the top,
      ! and edge 4 may not correspond to the left.
      read(IIN) numabsread,codeabsread(1),codeabsread(2),codeabsread(3), &
                codeabsread(4), typeabsread, ibegin_edge1(inum), iend_edge1(inum), &
                ibegin_edge2(inum), iend_edge2(inum), ibegin_edge3(inum), &
                iend_edge3(inum), ibegin_edge4(inum), iend_edge4(inum)

      if (numabsread < 1 .or. numabsread > nspec) &
        call exit_MPI(myrank,'Wrong absorbing element number')

      numabs(inum) = numabsread

      codeabs(IEDGE1,inum) = codeabsread(1) ! bottom
      codeabs(IEDGE2,inum) = codeabsread(2) ! right
      codeabs(IEDGE3,inum) = codeabsread(3) ! top
      codeabs(IEDGE4,inum) = codeabsread(4) ! left

      abs_boundary_type(inum) = typeabsread ! type == IBOTTOM / IRIGHT / ITOP / ILEFT

      ! check that a single edge is defined for each element cited
      ! (since elements with two absorbing edges MUST be cited twice, each time with a different "type" code
      if (count(codeabs(:,inum) .eqv. .true.) /= 1) then
        print *,'Error for absorbing element inum = ',inum
        call stop_the_code('must have one and only one absorbing edge per absorbing line cited')
      endif

    enddo

    ! detection of the corner element
    do inum = 1,num_abs_boundary_faces
      if (codeabs(IEDGE1,inum)) then
        ! bottom
        do inum_duplicate = 1,num_abs_boundary_faces
          if (inum == inum_duplicate) then
            ! left for blank, since no operation is needed.
            continue
          else
            if (numabs(inum) == numabs(inum_duplicate)) then
              if (codeabs(IEDGE4,inum_duplicate)) then
                ! left
                codeabs_corner(1,inum) = .true.
              endif
              if (codeabs(IEDGE2,inum_duplicate)) then
                ! right
                codeabs_corner(2,inum) = .true.
              endif
            endif
          endif
        enddo
      endif

      if (codeabs(IEDGE3,inum)) then
        ! top
        do inum_duplicate = 1,num_abs_boundary_faces
          if (inum == inum_duplicate) then
            ! left for blank, since no operation is needed.
            continue
          else
            if (numabs(inum) == numabs(inum_duplicate)) then
              if (codeabs(IEDGE4,inum_duplicate)) then
                ! left
                codeabs_corner(3,inum) = .true.
              endif
              if (codeabs(IEDGE2,inum_duplicate)) then
                ! right
                codeabs_corner(4,inum) = .true.
              endif
             endif
           endif
        enddo
      endif
    enddo

    ! detection of the corner element
    ! boundary element numbering
    do inum = 1,num_abs_boundary_faces
      if (codeabs(IEDGE1,inum)) then
        ! bottom
        nspec_bottom = nspec_bottom + 1
        ib_bottom(inum) =  nspec_bottom

      else if (codeabs(IEDGE2,inum)) then
        ! right
        nspec_right = nspec_right + 1
        ib_right(inum) =  nspec_right

      else if (codeabs(IEDGE3,inum)) then
        ! top
        nspec_top = nspec_top + 1
        ib_top(inum) = nspec_top

      else if (codeabs(IEDGE4,inum)) then
        ! left
        nspec_left = nspec_left + 1
        ib_left(inum) =  nspec_left

      else
        call stop_the_code('incorrect absorbing boundary element type read')
      endif
    enddo

  else
    ! if this MPI slice has no absorbing element at all
    num_abs_boundary_faces = 0
    nspec_left = 0
    nspec_right = 0
    nspec_bottom = 0
    nspec_top = 0
  endif

  ! sets up arrays for boundary routines
  if (num_abs_boundary_faces > 0) then
    allocate(abs_boundary_ispec(num_abs_boundary_faces),stat=ier)
    if (ier /= 0) stop 'error allocating array abs_boundary_ispec etc.'
  else
    ! dummy allocation
    allocate(abs_boundary_ispec(1))
  endif
  abs_boundary_ispec(:) = 0

  if (num_abs_boundary_faces > 0) then
    abs_boundary_ispec(:) = numabs(:)
  endif

  ! free memory
  deallocate(numabs)

  ! collects total values
  call sum_all_i(num_abs_boundary_faces, num_abs_boundary_faces_tot)
  call sum_all_i(nspec_left, nspec_left_tot)
  call sum_all_i(nspec_right, nspec_right_tot)
  call sum_all_i(nspec_bottom, nspec_bottom_tot)
  call sum_all_i(nspec_top, nspec_top_tot)

  ! user output
  if (PML_BOUNDARY_CONDITIONS .or. STACEY_ABSORBING_CONDITIONS) then
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Absorbing boundaries:'
      if (PML_BOUNDARY_CONDITIONS) &
        write(IMAIN,*) '  using PML boundary conditions'
      if (STACEY_ABSORBING_CONDITIONS) &
        write(IMAIN,*) '  using Stacey absorbing boundary conditions'
      ! for Stacey
      if (STACEY_ABSORBING_CONDITIONS) then
        write(IMAIN,*)
        write(IMAIN,*) '  Total number of absorbing elements: ',num_abs_boundary_faces_tot
        write(IMAIN,*) '    nspec_left   = ',nspec_left_tot
        write(IMAIN,*) '    nspec_right  = ',nspec_right_tot
        write(IMAIN,*) '    nspec_bottom = ',nspec_bottom_tot
        write(IMAIN,*) '    nspec_top    = ',nspec_top_tot
        write(IMAIN,*)
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  end subroutine read_mesh_databases_absorbing

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_acoustic_forcing()

! reads in absorbing edges

  use constants, only: IIN,IRIGHT,ILEFT,IBOTTOM,ITOP,IEDGE1,IEDGE2,IEDGE3,IEDGE4,IMAIN

  use specfem_par, only: myrank,nelem_acforcing,nspec,ACOUSTIC_FORCING, &
                         ibegin_edge1_acforcing,iend_edge1_acforcing,ibegin_edge2_acforcing,iend_edge2_acforcing, &
                         ibegin_edge3_acforcing,iend_edge3_acforcing,ibegin_edge4_acforcing,iend_edge4_acforcing, &
                         numacforcing,codeacforcing,typeacforcing, &
                         nspec_left_acforcing,nspec_right_acforcing,nspec_bottom_acforcing,nspec_top_acforcing, &
                         ib_right_acforcing,ib_left_acforcing,ib_bottom_acforcing,ib_top_acforcing

  implicit none

  ! local parameters
  integer :: inum,numacforcingread,typeacforcingread
  integer :: nelem_acforcing_all
  integer :: nspec_left_acforcing_all,nspec_right_acforcing_all,nspec_bottom_acforcing_all,nspec_top_acforcing_all
  integer :: ier
  logical :: codeacforcingread(4)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'reading acoustic forcing...'
    write(IMAIN,*) '  acoustic forcing                             = ',ACOUSTIC_FORCING
    write(IMAIN,*) '  number of acoustic forcing boundary elements = ',nelem_acforcing
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  if (ACOUSTIC_FORCING) then
    ! allocates arrays for acoustic forcing boundary conditions
    allocate(numacforcing(nelem_acforcing), &
             codeacforcing(4,nelem_acforcing), &
             typeacforcing(nelem_acforcing),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating acoustic forcing arrays')

    allocate(ibegin_edge1_acforcing(nelem_acforcing), &
             iend_edge1_acforcing(nelem_acforcing), &
             ibegin_edge3_acforcing(nelem_acforcing), &
             iend_edge3_acforcing(nelem_acforcing), &
             ibegin_edge4_acforcing(nelem_acforcing), &
             iend_edge4_acforcing(nelem_acforcing), &
             ibegin_edge2_acforcing(nelem_acforcing), &
             iend_edge2_acforcing(nelem_acforcing), &
             ib_left_acforcing(nelem_acforcing), &
             ib_right_acforcing(nelem_acforcing), &
             ib_bottom_acforcing(nelem_acforcing), &
             ib_top_acforcing(nelem_acforcing),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating acoustic forcing boundary arrays')
  else
    ! dummy allocation
    allocate(numacforcing(1),codeacforcing(1,1),typeacforcing(1))
    allocate(ibegin_edge1_acforcing(1),iend_edge1_acforcing(1),ibegin_edge3_acforcing(1),iend_edge3_acforcing(1))
    allocate(ibegin_edge4_acforcing(1),iend_edge4_acforcing(1),ibegin_edge2_acforcing(1),iend_edge2_acforcing(1))
    allocate(ib_left_acforcing(1),ib_right_acforcing(1),ib_bottom_acforcing(1),ib_top_acforcing(1))
  endif
  ! initializes
  codeacforcing(:,:) = .false.
  typeacforcing(:) = 0

  ibegin_edge1_acforcing(:) = 0
  iend_edge1_acforcing(:) = 0
  ibegin_edge3_acforcing(:) = 0
  iend_edge3_acforcing(:) = 0

  ibegin_edge4_acforcing(:) = 0
  iend_edge4_acforcing(:) = 0
  ibegin_edge2_acforcing(:) = 0
  iend_edge2_acforcing(:) = 0

  nspec_left_acforcing = 0
  nspec_right_acforcing = 0
  nspec_bottom_acforcing = 0
  nspec_top_acforcing = 0

  ib_right_acforcing(:) = 0
  ib_left_acforcing(:) = 0
  ib_bottom_acforcing(:) = 0
  ib_top_acforcing(:) = 0

  ! reads in forcing edges
  if (ACOUSTIC_FORCING) then

    ! reads forcing boundaries
    do inum = 1,nelem_acforcing

      ! beware here and below that external meshes (for instance coming from CUBIT or Gmsh)
      ! may have rotated elements and thus edge 1 may not correspond to the bottom,
      ! edge 2 may not correspond to the right, edge 3 may not correspond to the top,
      ! and edge 4 may not correspond to the left.
      read(IIN) numacforcingread,codeacforcingread(1),codeacforcingread(2),codeacforcingread(3), &
                codeacforcingread(4), typeacforcingread, ibegin_edge1_acforcing(inum), iend_edge1_acforcing(inum), &
                ibegin_edge2_acforcing(inum), iend_edge2_acforcing(inum), ibegin_edge3_acforcing(inum), &
                iend_edge3_acforcing(inum), ibegin_edge4_acforcing(inum), iend_edge4_acforcing(inum)

      ! checks index
      if (numacforcingread < 1 .or. numacforcingread > nspec) &
        call exit_MPI(myrank,'Wrong absorbing element number')

      numacforcing(inum) = numacforcingread

      codeacforcing(IEDGE1,inum) = codeacforcingread(1)
      codeacforcing(IEDGE2,inum) = codeacforcingread(2)
      codeacforcing(IEDGE3,inum) = codeacforcingread(3)
      codeacforcing(IEDGE4,inum) = codeacforcingread(4)

      typeacforcing(inum) = typeacforcingread

      ! check that a single edge is defined for each element cited
      ! (since elements with two absorbing edges MUST be cited twice, each time with a different "typeacforcing()" code
      if (count(codeacforcing(:,inum) .eqv. .true.) /= 1) then
        print *,'Error for absorbing element inum = ',inum
        call stop_the_code('must have one and only one absorbing edge per absorbing line cited')
      endif

    enddo

    ! boundary element numbering
    do inum = 1,nelem_acforcing
      if (typeacforcing(inum) == IBOTTOM) then
        nspec_bottom_acforcing = nspec_bottom_acforcing + 1
        ib_bottom_acforcing(inum) =  nspec_bottom_acforcing

      else if (typeacforcing(inum) == IRIGHT) then
        nspec_right_acforcing = nspec_right_acforcing + 1
        ib_right_acforcing(inum) =  nspec_right_acforcing

      else if (typeacforcing(inum) == ITOP) then
        nspec_top_acforcing = nspec_top_acforcing + 1
        ib_top_acforcing(inum) = nspec_top_acforcing

      else if (typeacforcing(inum) == ILEFT) then
        nspec_left_acforcing = nspec_left_acforcing + 1
        ib_left_acforcing(inum) =  nspec_left_acforcing

      else
        call stop_the_code('incorrect absorbing boundary element type read')
      endif
    enddo

    ! collects total number
    call sum_all_i(nelem_acforcing,nelem_acforcing_all)
    call sum_all_i(nspec_left_acforcing,nspec_left_acforcing_all)
    call sum_all_i(nspec_right_acforcing,nspec_right_acforcing_all)
    call sum_all_i(nspec_bottom_acforcing,nspec_bottom_acforcing_all)
    call sum_all_i(nspec_top_acforcing,nspec_top_acforcing_all)

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  Total number of acoustic forcing elements: ',nelem_acforcing_all
      write(IMAIN,*) '    nspec_left_acforcing   = ',nspec_left_acforcing_all
      write(IMAIN,*) '    nspec_right_acforcing  = ',nspec_right_acforcing_all
      write(IMAIN,*) '    nspec_bottom_acforcing = ',nspec_bottom_acforcing_all
      write(IMAIN,*) '    nspec_top_acforcing    = ',nspec_top_acforcing_all
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  end subroutine read_mesh_databases_acoustic_forcing

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_free_surf()

! reads acoustic free surface data

  use constants, only: IIN,IMAIN

  use specfem_par, only: myrank,nelem_acoustic_surface,acoustic_edges,acoustic_surface,any_acoustic_edges

  implicit none

  ! local parameters
  integer :: inum,nelem_acoustic_surface_all
  integer :: ier

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'reading free surface information...'
    write(IMAIN,*) '  number of acoustic free surface boundary elements = ',nelem_acoustic_surface
    call flush_IMAIN()
  endif

  ! sets acoustic edges flag
  if (nelem_acoustic_surface > 0) then
    any_acoustic_edges = .true.
  else
    any_acoustic_edges = .false.
  endif

  if (nelem_acoustic_surface > 0) then
    allocate(acoustic_edges(4,nelem_acoustic_surface), &
             acoustic_surface(5,nelem_acoustic_surface),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating acoustic free surface arrays')
  else
    ! dummy
    allocate(acoustic_edges(1,1),acoustic_surface(1,1),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating acoustic free surface arrays')
  endif
  ! initializes
  acoustic_edges(:,:) = 0

  ! reads in any possible free surface edges
  if (any_acoustic_edges) then
    do inum = 1,nelem_acoustic_surface
      read(IIN) acoustic_edges(1,inum), acoustic_edges(2,inum), acoustic_edges(3,inum), acoustic_edges(4,inum)
    enddo
  endif

  ! constructs (local) acoustic surface
  if (nelem_acoustic_surface > 0) then
    call construct_acoustic_surface ()
  endif

  ! collects total number
  call sum_all_i(nelem_acoustic_surface,nelem_acoustic_surface_all)

  if (nelem_acoustic_surface_all > 0) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  Total number of acoustic free surface elements: ',nelem_acoustic_surface_all
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  end subroutine read_mesh_databases_free_surf

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_coupled()

! reads acoustic elastic coupled edges
! reads acoustic poroelastic coupled edges
! reads poroelastic elastic coupled edges

  use constants, only: IIN,IMAIN

  use specfem_par, only: myrank

  use specfem_par, only: num_fluid_solid_edges,any_fluid_solid_edges, &
                         fluid_solid_acoustic_ispec,fluid_solid_elastic_ispec, &
                         fluid_solid_acoustic_iedge,fluid_solid_elastic_iedge, &
                         num_fluid_poro_edges,any_fluid_poro_edges, &
                         fluid_poro_acoustic_ispec,fluid_poro_poroelastic_ispec, &
                         fluid_poro_acoustic_iedge,fluid_poro_poroelastic_iedge, &
                         num_solid_poro_edges,any_solid_poro_edges, &
                         solid_poro_elastic_ispec,solid_poro_poroelastic_ispec, &
                         solid_poro_elastic_iedge,solid_poro_poroelastic_iedge

  implicit none

  ! local parameters
  integer :: inum
  integer :: fluid_solid_acoustic_ispec_read,fluid_solid_elastic_ispec_read, &
    fluid_poro_acoustic_ispec_read,fluid_poro_poro_ispec_read, &
    solid_poro_poro_ispec_read,solid_poro_elastic_ispec_read
  integer :: num_fluid_solid_edges_all,num_fluid_poro_edges_all,num_solid_poro_edges_all
  integer :: ier

  ! collects total number
  call sum_all_i(num_fluid_solid_edges,num_fluid_solid_edges_all)
  call sum_all_i(num_fluid_poro_edges,num_fluid_poro_edges_all)
  call sum_all_i(num_solid_poro_edges,num_solid_poro_edges_all)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'reading coupling surfaces...'
    write(IMAIN,*) '  number of fluid-solid edges  = ',num_fluid_solid_edges_all
    write(IMAIN,*) '  number of fluid-poro  edges  = ',num_fluid_poro_edges_all
    write(IMAIN,*) '  number of solid-poro  edges  = ',num_solid_poro_edges_all
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! sets flags fluid-solid domains
  if (num_fluid_solid_edges > 0) then
    any_fluid_solid_edges = .true.
  else
    any_fluid_solid_edges = .false.
  endif

  if (num_fluid_solid_edges > 0) then
    allocate(fluid_solid_acoustic_ispec(num_fluid_solid_edges), &
             fluid_solid_acoustic_iedge(num_fluid_solid_edges), &
             fluid_solid_elastic_ispec(num_fluid_solid_edges), &
             fluid_solid_elastic_iedge(num_fluid_solid_edges),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating fluid-solid arrays')
  else
    ! dummy allocation
    allocate(fluid_solid_acoustic_ispec(1),fluid_solid_acoustic_iedge(1), &
             fluid_solid_elastic_ispec(1),fluid_solid_elastic_iedge(1),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating fluid-solid arrays')
  endif
  ! initializes
  fluid_solid_acoustic_ispec(:) = 0; fluid_solid_acoustic_iedge(:) = 0
  fluid_solid_elastic_ispec(:) = 0; fluid_solid_elastic_iedge(:) = 0

  ! sets flags for poroelastic-acoustic coupled domains
  if (num_fluid_poro_edges > 0) then
    any_fluid_poro_edges = .true.
  else
    any_fluid_poro_edges = .false.
  endif

  if (num_fluid_poro_edges > 0) then
    allocate(fluid_poro_acoustic_ispec(num_fluid_poro_edges), &
             fluid_poro_acoustic_iedge(num_fluid_poro_edges), &
             fluid_poro_poroelastic_ispec(num_fluid_poro_edges), &
             fluid_poro_poroelastic_iedge(num_fluid_poro_edges),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating fluid-poro arrays')
  else
    ! dummy
    allocate(fluid_poro_acoustic_ispec(1),fluid_poro_acoustic_iedge(1), &
             fluid_poro_poroelastic_ispec(1),fluid_poro_poroelastic_iedge(1),stat=ier)
   if (ier /= 0) call stop_the_code('Error allocating fluid-poro arrays')
  endif
  ! initializes
  fluid_poro_acoustic_ispec(:) = 0; fluid_poro_acoustic_iedge(:) = 0
  fluid_poro_poroelastic_ispec(:) = 0; fluid_poro_poroelastic_iedge(:) = 0

  ! sets flags for poroelastic-solid coupled domains
  if (num_solid_poro_edges > 0) then
    any_solid_poro_edges = .true.
  else
    any_solid_poro_edges = .false.
  endif

  if (num_solid_poro_edges > 0) then
    allocate(solid_poro_elastic_ispec(num_solid_poro_edges), &
             solid_poro_elastic_iedge(num_solid_poro_edges), &
             solid_poro_poroelastic_ispec(num_solid_poro_edges), &
             solid_poro_poroelastic_iedge(num_solid_poro_edges),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating solid-poro arrays')
  else
    ! dummy
    allocate(solid_poro_elastic_ispec(1),solid_poro_elastic_iedge(1), &
             solid_poro_poroelastic_ispec(1),solid_poro_poroelastic_iedge(1),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating solid-poro arrays')
  endif
  ! initializes
  solid_poro_elastic_ispec(:) = 0; solid_poro_elastic_iedge(:) = 0
  solid_poro_poroelastic_ispec(:) = 0; solid_poro_poroelastic_iedge(:) = 0

  ! reads acoustic elastic coupled edges
  if (any_fluid_solid_edges) then
    do inum = 1, num_fluid_solid_edges
      read(IIN) fluid_solid_acoustic_ispec_read,fluid_solid_elastic_ispec_read

      fluid_solid_acoustic_ispec(inum) = fluid_solid_acoustic_ispec_read
      fluid_solid_elastic_ispec(inum) = fluid_solid_elastic_ispec_read
    enddo
  endif

  ! reads acoustic poroelastic coupled edges
  if (any_fluid_poro_edges) then
    do inum = 1, num_fluid_poro_edges
      read(IIN) fluid_poro_acoustic_ispec_read,fluid_poro_poro_ispec_read

      fluid_poro_acoustic_ispec(inum) = fluid_poro_acoustic_ispec_read
      fluid_poro_poroelastic_ispec(inum) = fluid_poro_poro_ispec_read
    enddo
  endif

  ! reads poroelastic elastic coupled edges
  if (any_solid_poro_edges) then
    do inum = 1, num_solid_poro_edges
      read(IIN) solid_poro_poro_ispec_read,solid_poro_elastic_ispec_read

      solid_poro_elastic_ispec(inum) = solid_poro_elastic_ispec_read
      solid_poro_poroelastic_ispec(inum) = solid_poro_poro_ispec_read
    enddo
  endif

  end subroutine read_mesh_databases_coupled

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_tangential()

! reads tangential detection curve

  use constants, only: IIN,IMAIN

  use specfem_par, only: myrank

  use specfem_par, only: nnodes_tangential_curve,nodes_tangential_curve, &
                         force_normal_to_surface,rec_normal_to_surface

  implicit none

  ! local parameters
  integer :: i,ier,nnodes_tangential_curve_all
  logical :: any_tangential_curve

  ! collects total number
  call sum_all_i(nnodes_tangential_curve,nnodes_tangential_curve_all)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'reading tangential curves...'
    write(IMAIN,*) '  number of tangential curve nodes = ',nnodes_tangential_curve_all
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! sets tangential flag
  if (nnodes_tangential_curve > 0) then
    any_tangential_curve = .true.
  else
    any_tangential_curve = .false.
  endif

  if (nnodes_tangential_curve > 0) then
    allocate(nodes_tangential_curve(2,nnodes_tangential_curve),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating tangential arrays')
  else
    ! dummy
    allocate(nodes_tangential_curve(1,1))
  endif
  ! initializes
  nodes_tangential_curve(:,:) = 0.d0

  ! reads tangential detection curve
  read(IIN) force_normal_to_surface,rec_normal_to_surface

  if (any_tangential_curve) then
    do i = 1, nnodes_tangential_curve
      read(IIN) nodes_tangential_curve(1,i),nodes_tangential_curve(2,i)
    enddo
  else
    force_normal_to_surface = .false.
    rec_normal_to_surface = .false.
  endif

  end subroutine read_mesh_databases_tangential

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_axial_elements()

! reads axial elements data

  use constants, only: IIN,IMAIN

  use specfem_par, only: myrank,nspec,nelem_on_the_axis,ispec_of_axial_elements,is_on_the_axis,AXISYM

  implicit none

  ! local parameters
  integer :: nelem_on_the_axis_total
  integer :: i,ier,ispec

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'reading axial elements...'
    write(IMAIN,*) '  number of axial elements     = ',nelem_on_the_axis
    call flush_IMAIN()
  endif

  ! axial flags
  allocate(is_on_the_axis(nspec),stat=ier)
  if (ier /= 0) call stop_the_code('Error: not enough memory to allocate array is_on_the_axis')

  is_on_the_axis(:) = .false.

  ! allocates ispec array for axial elements
  if (nelem_on_the_axis > 0) then
    allocate(ispec_of_axial_elements(nelem_on_the_axis),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating axial elements')
  else
    ! dummy array
    allocate(ispec_of_axial_elements(1))
  endif
  ! initializes
  ispec_of_axial_elements(:) = 0

  ! reads in any possible axial elements
  if (nelem_on_the_axis > 0) then
    do i = 1,nelem_on_the_axis
      read(IIN) ispec

      ! quick check
      if (ispec < 1 .or. ispec > nspec) call stop_the_code('Invalid ispec value, out of range in reading axial elements')

      ! stores element list
      ispec_of_axial_elements(i) = ispec
    enddo

    ! determines flags if element is on the axis
    call build_is_on_the_axis()
  endif

  ! collects total number of axial elements
  call sum_all_i(nelem_on_the_axis, nelem_on_the_axis_total)

  ! user output
  if (myrank == 0 .and. AXISYM) then
    write(IMAIN,*)
    write(IMAIN,*) ' Total number of elements on the axis : ',nelem_on_the_axis_total
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine read_mesh_databases_axial_elements

!
!-------------------------------------------------------------------------------------------------
!

  subroutine process_seismotype_line()

 ! This subroutine convert the string "seismotype" (e.g 2,3,6) to an integer array: "seismotypeVec"

  use specfem_par, only: NSIGTYPE,seismotype,seismotypeVec

  implicit none

    integer :: i
    character(len=512) :: stringCopy, stringCopy1
    character(len=52), parameter :: charVec = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ&
                                     &abcdefghijklmnopqrstuvwxyz'
    integer :: actual
    integer :: stringLen, integs
    integer :: res(512)

    interface
       subroutine AddToList(list, element)
         integer, dimension(:), allocatable, intent(inout) :: list
         integer, intent(in) :: element
       end subroutine AddToList
    end interface

    actual = 1
    NSIGTYPE = 0
    stringLen = len(charVec)
    stringCopy = trim(seismotype)

    ! Strip spaces :
    call StripSpaces(stringCopy)

    ! Strip chars :
    do while (actual < stringLen)
      call StripChar(stringCopy, charVec(actual:actual))
      actual = actual + 1
    enddo

    ! Convert to integer list :
    if (index(stringCopy, ',') < 0) then
      stop "This should never ever happen... this is an interesting bug really, memory broken somewhere?"
    else if (index(stringCopy, ',') == 0) then
      read(stringCopy, '(i3)') integs
      call AddToList(seismotypeVec, integs)
    else
      do while (index(stringCopy, ',') > 0)
        stringCopy1 = stringCopy(:index(stringCopy, ',')-1)
        read(stringCopy1, '(i3)') integs
        call AddToList(seismotypeVec, integs)
        stringCopy = stringCopy(index(stringCopy, ',')+1:)
        if (index(stringCopy, ',') < 1) then
          read(stringCopy, '(i3)') integs
          call AddToList(seismotypeVec, integs)
        endif
      enddo
    endif

    if (any(seismotypeVec == 0)) call stop_the_code("Seismotype line in Par_file contain a 0 or is malformed")

    ! Delete duplicates
    NSIGTYPE = 1
    res(:) = 0
    res(1) = seismotypeVec(1)
    do i = 2, size(seismotypeVec)
        ! if the number already exist in res check next
        if (any( res == seismotypeVec(i) )) cycle
        ! No match found so add it to the output
        NSIGTYPE = NSIGTYPE + 1
        res(NSIGTYPE) = seismotypeVec(i)
    enddo

    seismotypeVec = res(1:NSIGTYPE)

  end subroutine process_seismotype_line

!
!-------------------------------------------------------------------------------------------------
!

  subroutine StripChar(string,char)
    ! Remove all character "char" from string. Warning: "char" can't be " " !!!

    implicit none

    character(len=*),intent(inout) :: string
    character(len=512) :: stringCopy1, stringCopy2
    character,intent(in) :: char

    if (char == ' ') then
      stop 'This function can not be used to strip spaces, use StripSpaces instead'
    endif

    do while (index(string,char,back=.true.) > 0)
      stringCopy1 = string(:index(string,char,back=.true.)-1)
      stringCopy2 = string(index(string,char,back=.true.)+1:)
      string = trim(stringCopy1)//trim(stringCopy2)
    enddo

  end subroutine StripChar

!
!-------------------------------------------------------------------------------------------------
!

  subroutine StripSpaces(string)

  ! Remove all spaces from string

  implicit none

  character(len=*) :: string
  integer :: stringLen
  integer :: last, actual

  stringLen = len (string)
  last = 1
  actual = 1

  do while (actual < stringLen)
      if (string(last:last) == ' ') then
          actual = actual + 1
          string(last:last) = string(actual:actual)
          string(actual:actual) = ' '
      else
          last = last + 1
          if (actual < last) &
              actual = last
      endif
  enddo

  end subroutine StripSpaces

!
!-------------------------------------------------------------------------------------------------
!

  subroutine AddToList(list, element)

  ! Add an element to an integer list dynamically

  implicit none

  integer :: i, isize
  integer, intent(in) :: element
  integer, dimension(:), allocatable, intent(inout) :: list
  integer, dimension(:), allocatable :: clist


  if (allocated(list)) then
      isize = size(list)
      allocate(clist(isize+1))
      do i=1,isize
      clist(i) = list(i)
      enddo
      clist(isize+1) = element

      deallocate(list)
      call move_alloc(clist, list)

  else
      allocate(list(1))
      list(1) = element
  endif

end subroutine AddToList


