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

  use constants, only: IMAIN,IIN,DISPLAY_COLORS,DISPLAY_ELEMENT_NUMBERS_POSTSCRIPT,OUTPUT_FILES
  use specfem_par
  use specfem_par_noise, only: NOISE_TOMOGRAPHY
  use specfem_par_movie

  implicit none

  ! local parameters
  integer, external :: err_occurred
  integer :: ier,int_dummy
  character(len=MAX_STRING_LEN) :: prname, dummy

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
  open(unit=IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0 ) then
    if (myrank == 0) then
      print *,'Error opening database file: ',trim(prname)
      print *
      print *,'Please make sure that the mesher has been run before this solver simulation with the correct settings...'
    endif
    call exit_MPI(myrank,'Error opening file '//trim(OUTPUT_FILES)//'Database***.bin')
  endif

  !-------- starts reading init section

  !---  read job title and skip remaining titles of the input file
  read(IIN) simulation_title

  !---- read parameters from input file
  read(IIN) NOISE_TOMOGRAPHY, UNDO_ATTENUATION_AND_OR_PML

  read(IIN) nspec

  read(IIN) npgeo,nproc_read_from_database

  read(IIN) output_grid_Gnuplot,interpol

  read(IIN) NSTEP_BETWEEN_OUTPUT_INFO

  read(IIN) NSTEP_BETWEEN_OUTPUT_SEISMOS

  read(IIN) NSTEP_BETWEEN_OUTPUT_IMAGES

  read(IIN) PML_BOUNDARY_CONDITIONS

  read(IIN) ROTATE_PML_ACTIVATE

  read(IIN) ROTATE_PML_ANGLE

  read(IIN) K_MIN_PML

  read(IIN) K_MAX_PML

  read(IIN) damping_change_factor_acoustic

  read(IIN) damping_change_factor_elastic

  read(IIN) PML_PARAMETER_ADJUSTMENT

  read(IIN) read_external_mesh

  read(IIN) NELEM_PML_THICKNESS

  read(IIN) subsamp_seismos,imagetype_JPEG,imagetype_wavefield_dumps

  read(IIN) output_postscript_snapshot,output_color_image

  read(IIN) meshvect,modelvect,boundvect,cutsnaps,subsamp_postscript,sizemax_arrows

  read(IIN) anglerec

  read(IIN) initialfield

  read(IIN) add_Bielak_conditions_bottom,add_Bielak_conditions_right,add_Bielak_conditions_top,add_Bielak_conditions_left

  read(IIN) seismotype,imagetype_postscript

  read(IIN) MODEL

  read(IIN) SAVE_MODEL

  read(IIN) TOMOGRAPHY_FILE

  read(IIN) output_grid_ASCII,OUTPUT_ENERGY,NTSTEP_BETWEEN_OUTPUT_ENERGY,output_wavefield_dumps

  read(IIN) use_binary_for_wavefield_dumps

  read(IIN) ATTENUATION_VISCOELASTIC,ATTENUATION_PORO_FLUID_PART,ATTENUATION_VISCOACOUSTIC

  read(IIN) USE_SOLVOPT

  read(IIN) save_ASCII_seismograms

  read(IIN) save_binary_seismograms_single,save_binary_seismograms_double

  read(IIN) USE_TRICK_FOR_BETTER_PRESSURE

  read(IIN) COMPUTE_INTEGRATED_ENERGY_FIELD

  read(IIN) save_ASCII_kernels

  read(IIN) NSTEP_BETWEEN_COMPUTE_KERNELS

  read(IIN) NO_BACKWARD_RECONSTRUCTION

  read(IIN) DRAW_SOURCES_AND_RECEIVERS

  read(IIN) Q0_poroelastic,freq0_poroelastic

  read(IIN) AXISYM

  read(IIN) P_SV

  read(IIN) factor_subsample_image

  read(IIN) USE_CONSTANT_MAX_AMPLITUDE

  read(IIN) CONSTANT_MAX_AMPLITUDE_TO_USE

  read(IIN) USE_SNAPSHOT_NUMBER_IN_FILENAME

  read(IIN) DRAW_WATER_IN_BLUE

  read(IIN) US_LETTER

  read(IIN) POWER_DISPLAY_COLOR

  read(IIN) SU_FORMAT

  read(IIN) USER_T0

  read(IIN) time_stepping_scheme

  read(IIN) ADD_PERIODIC_CONDITIONS

  read(IIN) PERIODIC_HORIZ_DIST

  read(IIN) GPU_MODE

  read(IIN) setup_with_binary_database

  !---- read time step
  read(IIN) NSTEP,DT

  read(IIN) NT_DUMP_ATTENUATION

  ! read the ACOUSTIC_FORCING flag
  read(IIN) ACOUSTIC_FORCING

  ! 'NUMBER_OF_SIMULTANEOUS_RUNS'
  read(IIN) NUMBER_OF_SIMULTANEOUS_RUNS

  ! 'BROADCAST_SAME_MESH_AND_MODEL'
  read(IIN) BROADCAST_SAME_MESH_AND_MODEL

  !-------- finish reading init section
  ! sets time step for time scheme
  deltat = DT

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
    endif

    ! outputs parameters read
    write(IMAIN,200) npgeo,NDIM
    write(IMAIN,600) NSTEP_BETWEEN_OUTPUT_INFO,DISPLAY_COLORS,DISPLAY_ELEMENT_NUMBERS_POSTSCRIPT
    write(IMAIN,700) seismotype,anglerec
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
  'Display frequency . . . .(NSTEP_BETWEEN_OUTPUT_INFO) = ',i6/ 5x, &
  'Color display . . . . . . . . . . . . . . . (colors) = ',i6/ 5x, &
  ' == 0     black and white display              ',  / 5x, &
  ' == 1     color display                        ',  /5x, &
  'Numbered mesh . . . . . . . . . . . . . . .(numbers) = ',i6/ 5x, &
  ' == 0     do not number the mesh               ',  /5x, &
  ' == 1     number the mesh                      ')

700 format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Seismograms recording type . . . . . . .(seismotype) = ',i6/5x, &
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

#ifdef USE_MPI
  use mpi
#endif

  use constants, only: IIN,ADD_A_SMALL_CRACK_IN_THE_MEDIUM
  use specfem_par

  implicit none

  ! continues reading database files
  ! note: we opened the database file in read_mesh_for_init() and will continue reading from its current position
  !       thus, the ordering here must be consistent with the order in save_databases.f90

  ! reads in source infos
  call read_mesh_databases_sources()

  ! sets source parameters
  call set_source_parameters()

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

  end subroutine read_mesh_databases

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_sources()

! reads source parameters

  use constants, only: IIN
  use specfem_par
  use specfem_par_movie

  implicit none

  ! local parameters
  integer :: i_source,ier

  !----  read source information
  read(IIN) NSOURCES

  ! safety check
  ! note: in principle, the number of sources could be zero for noise simulations.
  !       however, we want to make sure to have one defined at least, even if not really needed.
  if (NSOURCES < 1) call stop_the_code('Need at least one source for running a simulation, please check...')

  ! allocates source information arrays
  allocate(source_type(NSOURCES), &
           time_function_type(NSOURCES), &
           name_of_source_file(NSOURCES), &
           burst_band_width(NSOURCES), &
           Mxx(NSOURCES), &
           Mxz(NSOURCES), &
           Mzz(NSOURCES), &
           f0_source(NSOURCES), &
           tshift_src(NSOURCES), &
           factor(NSOURCES), &
           anglesource(NSOURCES), &
           ispec_selected_source(NSOURCES), &
           iglob_source(NSOURCES), &
           source_courbe_eros(NSOURCES), &
           islice_selected_source(NSOURCES), &
           sourcearrays(NSOURCES,NDIM,NGLLX,NGLLZ),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating source arrays')

  ! source locations
  allocate(x_source(NSOURCES), &
           z_source(NSOURCES), &
           xi_source(NSOURCES), &
           gamma_source(NSOURCES),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating source arrays')

  ! initializes
  source_type(:) = 0
  time_function_type(:) = 0

  f0_source(:) = 0.d0
  tshift_src(:) = 0.d0
  factor(:) = 0.d0
  anglesource(:) = 0.d0

  Mxx(:) = 0.d0
  Mzz(:) = 0.d0
  Mxz(:) = 0.d0

  x_source(:) = 0.d0
  z_source(:) = 0.d0

  islice_selected_source(:) = 0
  ispec_selected_source(:) = 0
  iglob_source(:) = 0

  sourcearrays(:,:,:,:) = 0._CUSTOM_REAL

  ! reads in source info from Database file (check with routine save_databases_sources())
  do i_source = 1,NSOURCES
    read(IIN) source_type(i_source),time_function_type(i_source)
    read(IIN) name_of_source_file(i_source)
    read(IIN) burst_band_width(i_source)
    read(IIN) x_source(i_source),z_source(i_source)
    read(IIN) f0_source(i_source),tshift_src(i_source)
    read(IIN) factor(i_source),anglesource(i_source)
    read(IIN) Mxx(i_source),Mzz(i_source),Mxz(i_source)
  enddo

  !if (AXISYM) factor = factor/(TWO*PI)   !!!!! axisym TODO verify

  end subroutine read_mesh_databases_sources

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
  read(IIN) numat,ngnod,nspec,pointsdisp,plot_lowerleft_corner_only

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
    write(IMAIN,207) nspec_all,ngnod,NGLLX,NGLLZ,NGLLX*NGLLZ,pointsdisp,numat, &
                     nelem_acforcing_all,nelem_acoustic_surface_all
    call flush_IMAIN()
  endif

  ! allocates mesh arrays
  ! elements
  allocate(kmato(nspec))
  allocate(knods(ngnod,nspec))

  ! material
  allocate(density(2,numat))
  allocate(anisotropy(10,numat)) ! don't forget c22 value (it is used for AXISYM simulations only)

  allocate(porosity(numat), &
           tortuosity(numat), &
           permeability(3,numat),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating porosity,.. arrays')

  allocate(poroelastcoef(4,3,numat),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating poroelastcoef arrays')

  allocate(QKappa_attenuation(numat), &
           Qmu_attenuation(numat),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating attenuation arrays')


  ! output formats
107 format(/5x,'-- Spectral Elements --',//)

207 format(5x,'Number of spectral elements . . . . . . . . .  (nspec) =',i7,/5x, &
               'Number of control nodes per element . . . . . (ngnod) =',i7,/5x, &
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

  use constants, only: IIN,NGLLX,NGLLZ

  use specfem_par, only: N_SLS,ATTENUATION_f0_REFERENCE,READ_VELOCITIES_AT_f0

  implicit none

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

  use constants, only: IIN

  use specfem_par, only: nspec,ngnod,kmato,knods,region_CPML,f0_source

  implicit none

  ! local parameters
  integer :: n,k,ispec,kmato_read,pml_read
  integer :: ier
  integer, dimension(:), allocatable :: knods_read

  ! reads and sets the material properties
  call read_materials(f0_source(1))

  ! add support for using PML in MPI mode with external mesh
  allocate(region_CPML(nspec),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating region_CPML array')

  ! initializes
  kmato(:) = 0
  knods(:,:) = 0
  region_CPML(:) = 0

  ! temporary read array
  allocate(knods_read(ngnod),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating temporary knods array')

  ! reads spectral macrobloc data

  ! reads in values
  n = 0
  do ispec = 1,nspec
    ! format: #element_id  #material_id #node_id1 #node_id2 #...
    read(IIN) n,kmato_read,(knods_read(k),k = 1,ngnod),pml_read

    ! material association
    kmato(n) = kmato_read
    region_CPML(n) = pml_read

    ! element control node indices
    knods(:,n)= knods_read(:)
  enddo

  ! frees temporary array
  deallocate(knods_read)

  end subroutine read_mesh_databases_mato

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_interfaces()

! reads in interface dimensions

  use constants, only: IIN
  use specfem_par

  implicit none

  ! local parameters
  integer :: num_interface,ie,my_interfaces_read

  ! reads number of interfaces
  read(IIN) ninterface, max_interface_size

  ! allocates arrays for mpi/partition interfaces
  if (ninterface > 0) then
    allocate(my_neighbors(ninterface))
    allocate(my_nelmnts_neighbors(ninterface))
    allocate(my_interfaces(4,max_interface_size,ninterface))

    allocate(ibool_interfaces_acoustic(NGLLX*max_interface_size,ninterface))
    allocate(ibool_interfaces_elastic(NGLLX*max_interface_size,ninterface))
    allocate(ibool_interfaces_poroelastic(NGLLX*max_interface_size,ninterface))
    allocate(ibool_interfaces_ext_mesh_init(NGLLX*max_interface_size,ninterface))

    allocate(nibool_interfaces_acoustic(ninterface))
    allocate(nibool_interfaces_elastic(ninterface))
    allocate(nibool_interfaces_poroelastic(ninterface))
    allocate(nibool_interfaces_ext_mesh(ninterface))

    allocate(inum_interfaces_acoustic(ninterface))
    allocate(inum_interfaces_elastic(ninterface))
    allocate(inum_interfaces_poroelastic(ninterface))
  else
    allocate(my_neighbors(1))
    allocate(my_nelmnts_neighbors(1))
    allocate(my_interfaces(1,1,1))

    allocate(ibool_interfaces_acoustic(1,1))
    allocate(ibool_interfaces_elastic(1,1))
    allocate(ibool_interfaces_poroelastic(1,1))
    allocate(ibool_interfaces_ext_mesh_init(1,1))

    allocate(nibool_interfaces_acoustic(1))
    allocate(nibool_interfaces_elastic(1))
    allocate(nibool_interfaces_poroelastic(1))
    allocate(nibool_interfaces_ext_mesh(1))

    allocate(inum_interfaces_acoustic(1))
    allocate(inum_interfaces_elastic(1))
    allocate(inum_interfaces_poroelastic(1))
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

#ifdef USE_MPI
  use mpi
#endif

  use constants, only: IMAIN,IIN,IEDGE1,IEDGE2,IEDGE3,IEDGE4
  use specfem_par

  implicit none

  ! local parameters
  integer :: inum,inum_duplicate,numabsread,typeabsread
  logical :: codeabsread(4)

  integer :: nelemabs_tot,nspec_left_tot,nspec_right_tot,nspec_bottom_tot,nspec_top_tot
  integer :: ier

  ! determines flag for absorbing boundaries
  if (nelemabs <= 0) then
    anyabs = .false.
    ! uses dummy value for allocation
    nelemabs = 1
  else
    anyabs = .true.
  endif

  ! sets Stacey flag
  if (anyabs .and. (.not. PML_BOUNDARY_CONDITIONS)) then
    STACEY_ABSORBING_CONDITIONS = .true.
  else
    STACEY_ABSORBING_CONDITIONS = .false.
  endif

  ! allocate arrays for absorbing boundary conditions
  allocate(numabs(nelemabs), &
           codeabs(4,nelemabs),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating absorbing arrays')

  !---codeabs_corner(1,nelemabs) denotes whether element is on bottom-left corner of absorbing boundary or not
  !---codeabs_corner(2,nelemabs) denotes whether element is on bottom-right corner of absorbing boundary or not
  !---codeabs_corner(3,nelemabs) denotes whether element is on top-left corner of absorbing boundary or not
  !---codeabs_corner(4,nelemabs) denotes whether element is on top-right corner of absorbing boundary or not
  allocate(codeabs_corner(4,nelemabs),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating absorbing codeabs_corner array')

  allocate(typeabs(nelemabs))

  allocate(ibegin_edge1(nelemabs))
  allocate(iend_edge1(nelemabs))

  allocate(ibegin_edge2(nelemabs))
  allocate(iend_edge2(nelemabs))

  allocate(ibegin_edge3(nelemabs))
  allocate(iend_edge3(nelemabs))

  allocate(ibegin_edge4(nelemabs))
  allocate(iend_edge4(nelemabs))

  ! poroelastic
  allocate(ibegin_edge1_poro(nelemabs))
  allocate(iend_edge1_poro(nelemabs))

  allocate(ibegin_edge2_poro(nelemabs))
  allocate(iend_edge2_poro(nelemabs))

  allocate(ibegin_edge3_poro(nelemabs))
  allocate(iend_edge3_poro(nelemabs))

  allocate(ibegin_edge4_poro(nelemabs))
  allocate(iend_edge4_poro(nelemabs))

  allocate(ib_left(nelemabs), &
           ib_right(nelemabs), &
           ib_bottom(nelemabs), &
           ib_top(nelemabs),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating absorbing boundary arrays')

  ! initializes
  codeabs(:,:) = .false.
  codeabs_corner(:,:) = .false.
  typeabs(:) = 0

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
    do inum = 1,nelemabs

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

      codeabs(IEDGE1,inum) = codeabsread(1)
      codeabs(IEDGE2,inum) = codeabsread(2)
      codeabs(IEDGE3,inum) = codeabsread(3)
      codeabs(IEDGE4,inum) = codeabsread(4)

      typeabs(inum) = typeabsread

      ! check that a single edge is defined for each element cited
      ! (since elements with two absorbing edges MUST be cited twice, each time with a different "typeabs()" code
      if (count(codeabs(:,inum) .eqv. .true.) /= 1) then
        print *,'Error for absorbing element inum = ',inum
        call stop_the_code('must have one and only one absorbing edge per absorbing line cited')
      endif

    enddo

    ! detection of the corner element
    do inum = 1,nelemabs
      if (codeabs(IEDGE1,inum)) then
        ! bottom
        do inum_duplicate = 1,nelemabs
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
        do inum_duplicate = 1,nelemabs
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
    do inum = 1,nelemabs
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
    nelemabs = 0
    nspec_left = 0
    nspec_right = 0
    nspec_bottom = 0
    nspec_top = 0
  endif

  ! collects total values
  call sum_all_i(nelemabs, nelemabs_tot)
  call sum_all_i(nspec_left, nspec_left_tot)
  call sum_all_i(nspec_right, nspec_right_tot)
  call sum_all_i(nspec_bottom, nspec_bottom_tot)
  call sum_all_i(nspec_top, nspec_top_tot)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    if (PML_BOUNDARY_CONDITIONS .or. STACEY_ABSORBING_CONDITIONS) then
      write(IMAIN,*) 'Absorbing boundaries:'
      if (PML_BOUNDARY_CONDITIONS) &
        write(IMAIN,*) '  using PML boundary conditions'
      if (STACEY_ABSORBING_CONDITIONS) &
        write(IMAIN,*) '  using Stacey absorbing boundary conditions'
      ! for Stacey
      if (STACEY_ABSORBING_CONDITIONS) then
        write(IMAIN,*)
        write(IMAIN,*) 'Number of absorbing elements: ',nelemabs_tot
        write(IMAIN,*) '  nspec_left   = ',nspec_left_tot
        write(IMAIN,*) '  nspec_right  = ',nspec_right_tot
        write(IMAIN,*) '  nspec_bottom = ',nspec_bottom_tot
        write(IMAIN,*) '  nspec_top    = ',nspec_top_tot
        write(IMAIN,*)
      endif
    endif
    call flush_IMAIN()
  endif

  ! allocates arrays
  if (anyabs) then
    ! files to save absorbed waves needed to reconstruct backward wavefield for adjoint method
    if (any_elastic .and. (SAVE_FORWARD .or. SIMULATION_TYPE == 3) .and. STACEY_ABSORBING_CONDITIONS) then
      allocate(b_absorb_elastic_left(NDIM,NGLLZ,nspec_left,NSTEP))
      allocate(b_absorb_elastic_right(NDIM,NGLLZ,nspec_right,NSTEP))
      allocate(b_absorb_elastic_bottom(NDIM,NGLLX,nspec_bottom,NSTEP))
      allocate(b_absorb_elastic_top(NDIM,NGLLX,nspec_top,NSTEP))
    else
      allocate(b_absorb_elastic_left(1,1,1,1))
      allocate(b_absorb_elastic_right(1,1,1,1))
      allocate(b_absorb_elastic_bottom(1,1,1,1))
      allocate(b_absorb_elastic_top(1,1,1,1))
    endif
    if (any_poroelastic .and. (SAVE_FORWARD .or. SIMULATION_TYPE == 3) .and. STACEY_ABSORBING_CONDITIONS) then
      allocate(b_absorb_poro_s_left(NDIM,NGLLZ,nspec_left,NSTEP))
      allocate(b_absorb_poro_s_right(NDIM,NGLLZ,nspec_right,NSTEP))
      allocate(b_absorb_poro_s_bottom(NDIM,NGLLX,nspec_bottom,NSTEP))
      allocate(b_absorb_poro_s_top(NDIM,NGLLX,nspec_top,NSTEP))
      allocate(b_absorb_poro_w_left(NDIM,NGLLZ,nspec_left,NSTEP))
      allocate(b_absorb_poro_w_right(NDIM,NGLLZ,nspec_right,NSTEP))
      allocate(b_absorb_poro_w_bottom(NDIM,NGLLX,nspec_bottom,NSTEP))
      allocate(b_absorb_poro_w_top(NDIM,NGLLX,nspec_top,NSTEP))
    else
      allocate(b_absorb_poro_s_left(1,1,1,1))
      allocate(b_absorb_poro_s_right(1,1,1,1))
      allocate(b_absorb_poro_s_bottom(1,1,1,1))
      allocate(b_absorb_poro_s_top(1,1,1,1))
      allocate(b_absorb_poro_w_left(1,1,1,1))
      allocate(b_absorb_poro_w_right(1,1,1,1))
      allocate(b_absorb_poro_w_bottom(1,1,1,1))
      allocate(b_absorb_poro_w_top(1,1,1,1))
    endif
    if (any_acoustic .and. (SAVE_FORWARD .or. SIMULATION_TYPE == 3) .and. STACEY_ABSORBING_CONDITIONS) then
      allocate(b_absorb_acoustic_left(NGLLZ,nspec_left,NSTEP))
      allocate(b_absorb_acoustic_right(NGLLZ,nspec_right,NSTEP))
      allocate(b_absorb_acoustic_bottom(NGLLX,nspec_bottom,NSTEP))
      allocate(b_absorb_acoustic_top(NGLLX,nspec_top,NSTEP))
    else
      allocate(b_absorb_acoustic_left(1,1,1))
      allocate(b_absorb_acoustic_right(1,1,1))
      allocate(b_absorb_acoustic_bottom(1,1,1))
      allocate(b_absorb_acoustic_top(1,1,1))
    endif

  else
    ! dummy arrays
    ! elastic domains
    if (.not. allocated(b_absorb_elastic_left)) then
      allocate(b_absorb_elastic_left(1,1,1,1))
      allocate(b_absorb_elastic_right(1,1,1,1))
      allocate(b_absorb_elastic_bottom(1,1,1,1))
      allocate(b_absorb_elastic_top(1,1,1,1))
    endif
    ! poroelastic domains
    if (.not. allocated(b_absorb_poro_s_left)) then
      allocate(b_absorb_poro_s_left(1,1,1,1))
      allocate(b_absorb_poro_s_right(1,1,1,1))
      allocate(b_absorb_poro_s_bottom(1,1,1,1))
      allocate(b_absorb_poro_s_top(1,1,1,1))
      allocate(b_absorb_poro_w_left(1,1,1,1))
      allocate(b_absorb_poro_w_right(1,1,1,1))
      allocate(b_absorb_poro_w_bottom(1,1,1,1))
      allocate(b_absorb_poro_w_top(1,1,1,1))
    endif
    ! acoustic domains
    if (.not. allocated(b_absorb_acoustic_left)) then
      allocate(b_absorb_acoustic_left(1,1,1))
      allocate(b_absorb_acoustic_right(1,1,1))
      allocate(b_absorb_acoustic_bottom(1,1,1))
      allocate(b_absorb_acoustic_top(1,1,1))
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

  if (.not. ACOUSTIC_FORCING) then
    ! dummy value for allocation
    nelem_acforcing = 1
  endif

  ! allocates arrays for acoustic forcing boundary conditions
  allocate(numacforcing(nelem_acforcing), &
           codeacforcing(4,nelem_acforcing), &
           typeacforcing(nelem_acforcing),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating acoustic forcing arrays')

  allocate(ibegin_edge1_acforcing(nelem_acforcing))
  allocate(iend_edge1_acforcing(nelem_acforcing))
  allocate(ibegin_edge3_acforcing(nelem_acforcing))
  allocate(iend_edge3_acforcing(nelem_acforcing))

  allocate(ibegin_edge4_acforcing(nelem_acforcing))
  allocate(iend_edge4_acforcing(nelem_acforcing))
  allocate(ibegin_edge2_acforcing(nelem_acforcing))
  allocate(iend_edge2_acforcing(nelem_acforcing))

  allocate(ib_left_acforcing(nelem_acforcing), &
           ib_right_acforcing(nelem_acforcing), &
           ib_bottom_acforcing(nelem_acforcing), &
           ib_top_acforcing(nelem_acforcing),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating acoustic forcing boundary arrays')

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
      write(IMAIN,*) 'Number of acoustic forcing elements: ',nelem_acforcing_all
      write(IMAIN,*) '  nspec_left_acforcing = ',nspec_left_acforcing_all
      write(IMAIN,*) '  nspec_right_acforcing = ',nspec_right_acforcing_all
      write(IMAIN,*) '  nspec_bottom_acforcing = ',nspec_bottom_acforcing_all
      write(IMAIN,*) '  nspec_top_acforcing = ',nspec_top_acforcing_all
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

  ! sets acoustic edges flag
  if (nelem_acoustic_surface > 0) then
    any_acoustic_edges = .true.
  else
    any_acoustic_edges = .false.
    ! dummy value for allocation
    nelem_acoustic_surface = 1
  endif

  allocate(acoustic_edges(4,nelem_acoustic_surface), &
           acoustic_surface(5,nelem_acoustic_surface),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating acoustic free surface arrays')

  ! initializes
  acoustic_edges(:,:) = 0

  ! reads in any possible free surface edges
  if (any_acoustic_edges) then
    do inum = 1,nelem_acoustic_surface
      read(IIN) acoustic_edges(1,inum), acoustic_edges(2,inum), acoustic_edges(3,inum), acoustic_edges(4,inum)
    enddo
  endif

  ! resets nelem_acoustic_surface
  if (any_acoustic_edges .eqv. .false. ) nelem_acoustic_surface = 0

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
      write(IMAIN,*) 'Number of free surface elements: ',nelem_acoustic_surface_all
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

  use constants, only: IIN

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
  integer :: ier

  ! sets flags fluid-solid domains
  if (num_fluid_solid_edges > 0) then
    any_fluid_solid_edges = .true.
  else
    any_fluid_solid_edges = .false.
    ! dummy value for allocation
    num_fluid_solid_edges = 1
  endif

  allocate(fluid_solid_acoustic_ispec(num_fluid_solid_edges), &
           fluid_solid_acoustic_iedge(num_fluid_solid_edges), &
           fluid_solid_elastic_ispec(num_fluid_solid_edges), &
           fluid_solid_elastic_iedge(num_fluid_solid_edges),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating fluid-solid arrays')

  ! sets flags for poroelastic-acoustic coupled domains
  if (num_fluid_poro_edges > 0) then
    any_fluid_poro_edges = .true.
  else
    any_fluid_poro_edges = .false.
    num_fluid_poro_edges = 1
  endif

  allocate(fluid_poro_acoustic_ispec(num_fluid_poro_edges))
  allocate(fluid_poro_acoustic_iedge(num_fluid_poro_edges))
  allocate(fluid_poro_poroelastic_ispec(num_fluid_poro_edges))
  allocate(fluid_poro_poroelastic_iedge(num_fluid_poro_edges))

  ! sets flags for poroelastic-solid coupled domains
  if (num_solid_poro_edges > 0) then
    any_solid_poro_edges = .true.
  else
    any_solid_poro_edges = .false.
    num_solid_poro_edges = 1
  endif

  allocate(solid_poro_elastic_ispec(num_solid_poro_edges))
  allocate(solid_poro_elastic_iedge(num_solid_poro_edges))
  allocate(solid_poro_poroelastic_ispec(num_solid_poro_edges))
  allocate(solid_poro_poroelastic_iedge(num_solid_poro_edges))

  ! initializes
  fluid_solid_acoustic_ispec(:) = 0
  fluid_solid_elastic_ispec(:) = 0

  fluid_poro_acoustic_ispec(:) = 0
  fluid_poro_poroelastic_ispec(:) = 0
  solid_poro_elastic_ispec(:) = 0
  solid_poro_poroelastic_ispec(:) = 0

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

  ! resets counters
  if (any_fluid_solid_edges .eqv. .false. ) num_fluid_solid_edges = 0
  if (any_fluid_poro_edges .eqv. .false. ) num_fluid_poro_edges = 0
  if (any_solid_poro_edges .eqv. .false. ) num_solid_poro_edges = 0

  end subroutine read_mesh_databases_coupled

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_tangential()

! reads tangential detection curve

  use constants, only: IIN

  use specfem_par, only: nnodes_tangential_curve,nodes_tangential_curve, &
                          force_normal_to_surface,rec_normal_to_surface, &
                          any_tangential_curve,dist_tangential_detection_curve

  implicit none

  ! local parameters
  integer :: i,ier

  ! sets tangential flag
  if (nnodes_tangential_curve > 0) then
    any_tangential_curve = .true.
  else
    any_tangential_curve = .false.
    ! dummy value for allocation
    nnodes_tangential_curve = 1
  endif

  allocate(nodes_tangential_curve(2,nnodes_tangential_curve), &
           dist_tangential_detection_curve(nnodes_tangential_curve),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating tangential arrays')

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

  ! resets nnode_tangential_curve
  if (any_tangential_curve .eqv. .false. ) nnodes_tangential_curve = 0

  end subroutine read_mesh_databases_tangential

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_axial_elements()

! reads axial elements data

#ifdef USE_MPI
  use mpi
#endif

  use constants, only: IIN,IMAIN

  use specfem_par, only: myrank,nspec,nelem_on_the_axis,ispec_of_axial_elements,is_on_the_axis,AXISYM

  implicit none

  ! local parameters
  integer :: nelem_on_the_axis_total
  integer :: i,ier,ispec

  ! axial flags
  allocate(is_on_the_axis(nspec),stat=ier)
  if (ier /= 0) call stop_the_code('Error: not enough memory to allocate array is_on_the_axis')

  is_on_the_axis(:) = .false.

  ! allocates ispec array for axial elements
  if (nelem_on_the_axis == 0) then
    ! dummy array
    allocate(ispec_of_axial_elements(1))
  else
    allocate(ispec_of_axial_elements(nelem_on_the_axis))
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
    write(IMAIN,*) 'Total number of elements on the axis : ',nelem_on_the_axis_total
    call flush_IMAIN()
  endif

  end subroutine read_mesh_databases_axial_elements

