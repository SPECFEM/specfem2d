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

  module parameter_file_par

  ! note: we use this module definition only to be able to allocate
  !          arrays for receiverlines and materials in this subroutine rather than in the main
  !          routine in meshfem2D.F90

  ! note 2: the filename ending is .F90 to have pre-compilation with pragmas
  !            (like #ifndef USE_MPI) working properly

  implicit none

  !--------------------------------------------------------------
  ! variables for preparation for compuation (simulation type, mesher, et al)
  !--------------------------------------------------------------
  ! name assigned to the running example
  character(len=100) :: title

  ! simulation type
  integer :: SIMULATION_TYPE,NOISE_TOMOGRAPHY
  logical :: SAVE_FORWARD,UNDO_ATTENUATION
  logical :: AXISYM
  logical :: P_SV

  ! computational platform type
  logical :: GPU_MODE

  ! input file name of TOMOGRAPHY
  character(len=100) :: TOMOGRAPHY_FILE

  ! input parameter for inner mesher
  double precision :: xmin_param,xmax_param
  integer :: nx_param
  integer :: ngnod
  character(len=100) :: interfacesfile
  logical, dimension(:),allocatable :: record_at_surface_same_vertical

  ! mesh files when using external mesh
  character(len=100) :: MODEL, SAVE_MODEL

  logical :: read_external_mesh

  character(len=256) :: mesh_file, nodes_coords_file, materials_file, &
                        free_surface_file, acoustic_forcing_surface_file, &
                        absorbing_surface_file, CPML_element_file, &
                        axial_elements_file

  character(len=256)  :: tangential_detection_curve_file

  ! material file for changing the model parameter for inner mesh or updating the
  ! the material for an existed mesh
  logical :: assign_external_model, READ_EXTERNAL_SEP_FILE

  !--------------------------------------------------------------
  ! variables for compuation
  !--------------------------------------------------------------
  ! variables used for partitioning
  integer :: NPROC, partitioning_method

  ! variables used for attenuation
  integer :: N_SLS
  logical :: READ_VELOCITIES_AT_f0,ATTENUATION_VISCOELASTIC_SOLID,ATTENUATION_PORO_FLUID_PART
  double precision  :: f0_attenuation
  double precision :: Q0,freq0

  ! variables used for plane wave incidence
  logical :: initialfield
  logical :: add_Bielak_conditions_left,add_Bielak_conditions_right,add_Bielak_conditions_bottom,add_Bielak_conditions_top, &
             add_Bielak_conditions

  ! variables used for absorbing boundary condition
  logical :: STACEY_ABSORBING_CONDITIONS,absorbbottom,absorbright,absorbtop,absorbleft

  integer :: NELEM_PML_THICKNESS
  logical :: PML_BOUNDARY_CONDITIONS,ROTATE_PML_ACTIVATE
  double precision :: ROTATE_PML_ANGLE

  ! for horizontal periodic conditions: detect common points between left and right edges
  logical :: ADD_PERIODIC_CONDITIONS
  ! horizontal periodicity distance for periodic conditions
  double precision :: PERIODIC_HORIZ_DIST

  ! variables used for source-receiver geometry
  integer :: NSOURCES

  ! acoustic forcing of an acoustic medium at a rigid interface
  logical :: ACOUSTIC_FORCING

  logical :: force_normal_to_surface
  logical :: use_existing_STATIONS

  integer :: nreceiversets
  double precision :: anglerec
  logical :: rec_normal_to_surface
  integer, dimension(:),allocatable :: nrec
  double precision, dimension(:),allocatable :: xdeb,zdeb,xfin,zfin

  ! variables used for iteration
  integer :: NT_DUMP_ATTENUATION

  ! time steps
  integer :: NSTEP
  double precision :: DT

  ! value of time_stepping_scheme to decide which time scheme will be used
  ! # 1 = Newmark (2nd order), 2 = LDDRK4-6 (4th-order 6-stage low storage Runge-Kutta)
  ! 3 = classical 4th-order 4-stage Runge-Kutta
  integer :: time_stepping_scheme

  ! use this t0 as earliest starting time rather than the automatically calculated one
  ! (must be positive and bigger than the automatically one to be effective;
  !  simulation will start at t = - t0)
  double precision :: USER_T0

  ! to store density and velocity model
  integer :: nbmodels
  integer, dimension(:),allocatable :: num_material

  integer, dimension(:),allocatable :: icodemat

  double precision, dimension(:),allocatable :: rho_s,cp,cs, &
    aniso3,aniso4,aniso5,aniso6,aniso7,aniso8,aniso9,aniso10,aniso11,aniso12,QKappa,Qmu

  double precision, dimension(:),allocatable :: rho_f,phi,tortuosity,permxx,permxz,&
       permzz,kappa_s,kappa_f,kappa_fr,eta_f,mu_fr

  !--------------------------------------------------------------
  ! variables used for output
  !--------------------------------------------------------------
  ! for later check of the grid
  logical :: output_grid_Gnuplot,output_grid_ASCII

  ! general information during the computation
  integer :: NSTEP_BETWEEN_OUTPUT_INFO

  ! for plotting the curve of energy
  logical :: output_energy

  ! kernel output in case of adjoint simulation
  logical :: save_ASCII_kernels

  ! seismogram
  integer :: seismotype,NSTEP_BETWEEN_OUTPUT_SEISMOS
  integer :: subsamp_seismos
  logical :: save_ASCII_seismograms,save_binary_seismograms_single,save_binary_seismograms_double
  logical :: USE_TRICK_FOR_BETTER_PRESSURE
  ! output seismograms in Seismic Unix format (adjoint traces will be read in the same format)
  logical :: SU_FORMAT

  ! Integrated energy field output
  logical :: COMPUTE_INTEGRATED_ENERGY_FIELD
  ! wave field
  integer :: NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS
  logical :: output_wavefield_dumps,use_binary_for_wavefield_dumps
  integer :: imagetype_wavefield_dumps

  ! image
  logical :: output_color_image
  integer :: imagetype_JPEG
  logical :: DRAW_SOURCES_AND_RECEIVERS
  ! use snapshot number in the file name of JPG color snapshots instead of the time step
  logical :: USE_SNAPSHOT_NUMBER_IN_FILENAME
  ! factor to subsample color images output by the code (useful for very large models)
  integer :: NSTEP_BETWEEN_OUTPUT_IMAGES
  double precision :: factor_subsample_image
  ! by default the code normalizes each image independently to its maximum; use this option to use the global maximum below instead
  logical :: USE_CONSTANT_MAX_AMPLITUDE
  ! constant maximum amplitude to use for all color images if the USE_CONSTANT_MAX_AMPLITUDE option is true
  double precision :: CONSTANT_MAX_AMPLITUDE_TO_USE
  logical :: plot_lowerleft_corner_only
  ! display acoustic layers as constant blue, because they likely correspond to water in the case of ocean acoustics
  ! or in the case of offshore oil industry experiments.
  ! (if off, display them as greyscale, as for elastic or poroelastic elements)
  logical :: DRAW_WATER_IN_BLUE
  ! non linear display to enhance small amplitudes in color images
  double precision :: POWER_DISPLAY_COLOR

  ! postscript
  logical :: output_postscript_snapshot
  integer :: imagetype_postscript
  double precision :: cutsnaps
  logical :: meshvect,modelvect,boundvect,interpol
  integer :: pointsdisp,subsamp_postscript
  double precision :: sizemax_arrows

  ! US letter paper or European A4
  logical :: US_LETTER

  end module parameter_file_par

!
!---------------------------------------------------------------------------------------
!

  module source_file_par

  implicit none

  ! source type parameters
  integer, dimension(:),pointer ::  source_type,time_function_type
  ! location
  double precision, dimension(:),pointer :: xs,zs
  ! moment tensor
  double precision, dimension(:),pointer :: Mxx,Mzz,Mxz
  ! source parameters
  double precision, dimension(:),pointer :: f0_source,tshift_src,anglesource,factor,burst_band_width
  ! flag for fixation to surface
  logical, dimension(:),allocatable ::  source_surf
  ! File name can't exceed 100 characters
  character(len=100), dimension(:),allocatable :: name_of_source_file

  end module source_file_par

!
!---------------------------------------------------------------------------------------
!

  module decompose_par

  implicit none

  ! variables used for storing info about the mesh and partitions
  integer, dimension(:), allocatable  :: my_interfaces
  integer, dimension(:), allocatable  :: my_nb_interfaces

  end module decompose_par

!
!---------------------------------------------------------------------------------------
!

  module part_unstruct_par

! This module contains subroutines related to unstructured meshes and partitioning of the
! corresponding graphs.

  implicit none

  integer :: nelmnts
  integer, dimension(:), allocatable  :: elmnts
  integer, dimension(:), allocatable  :: elmnts_bis
  integer, dimension(:), allocatable  :: vwgt
  integer, dimension(:), allocatable  :: glob2loc_elmnts
  integer, dimension(:), allocatable  :: part

  integer :: nb_edges
  integer, dimension(:), allocatable  :: adjwgt

  integer, dimension(:), allocatable  :: xadj_g
  integer, dimension(:), allocatable  :: adjncy_g

  integer :: nnodes
  double precision, dimension(:,:), allocatable  :: nodes_coords
  integer, dimension(:), allocatable  :: nnodes_elmnts
  integer, dimension(:), allocatable  :: nodes_elmnts
  integer, dimension(:), allocatable  :: glob2loc_nodes_nparts
  integer, dimension(:), allocatable  :: glob2loc_nodes_parts
  integer, dimension(:), allocatable  :: glob2loc_nodes

  ! interface data
  integer :: ninterfaces
  integer, dimension(:), allocatable  :: tab_size_interfaces, tab_interfaces

  integer :: nelem_acoustic_surface
  integer, dimension(:,:), allocatable  :: acoustic_surface
  integer :: nelem_acoustic_surface_loc

  integer :: nelem_on_the_axis
  integer, dimension(:), allocatable  :: ispec_of_axial_elements
  integer, dimension(:), allocatable  :: inode1_axial_elements, inode2_axial_elements
  integer :: nelem_on_the_axis_loc

  integer :: nelemabs
  integer, dimension(:,:), allocatable  :: abs_surface
  logical, dimension(:,:), allocatable  :: abs_surface_char
  integer, dimension(:), allocatable  :: abs_surface_merge,abs_surface_type
  integer :: nelemabs_loc

  integer :: nelemabs_merge
  integer, dimension(:), allocatable  :: ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
       ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2

  ! for acoustic/elastic coupled elements
  integer :: nedges_coupled
  integer, dimension(:,:), allocatable  :: edges_coupled

  ! for acoustic/poroelastic coupled elements
  integer :: nedges_acporo_coupled
  integer, dimension(:,:), allocatable  :: edges_acporo_coupled

  ! for poroelastic/elastic coupled elements
  integer :: nedges_elporo_coupled
  integer, dimension(:,:), allocatable  :: edges_elporo_coupled

  ! for acoustic forcing elements
  integer :: nelemacforcing
  integer, dimension(:,:), allocatable :: acforcing_surface
  logical, dimension(:,:), allocatable  :: acforcing_surface_char
  integer, dimension(:), allocatable  :: acforcing_surface_merge,acforcing_surface_type
  integer :: nelemacforcing_loc

  integer :: nelemacforcing_merge
  integer, dimension(:), allocatable  :: ibegin_edge1_acforcing,iend_edge1_acforcing, &
       ibegin_edge3_acforcing,iend_edge3_acforcing,ibegin_edge4_acforcing,iend_edge4_acforcing, &
       ibegin_edge2_acforcing,iend_edge2_acforcing

  ! vertical layers
  integer :: number_of_layers
  integer, dimension(:), allocatable :: nz_layer

  ! variables used for tangential detection
  integer ::  nnodes_tangential_curve
  double precision, dimension(:,:), allocatable  :: nodes_tangential_curve

  ! interface file data
  integer :: max_npoints_interface,number_of_interfaces
  integer :: nx,nz,nxread,nzread

  ! coordinates of the grid points of the mesh
  double precision, dimension(:,:), allocatable :: grid_point_x,grid_point_z

  integer :: npoints_interface_top
  double precision, dimension(:), allocatable :: xinterface_top,zinterface_top,coefs_interface_top

  end module part_unstruct_par
