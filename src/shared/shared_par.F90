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

! note: the filename ending is .F90 to have pre-compilation with pragmas
!            (like #ifndef USE_MPI) working properly

module constants

  include "constants.h"

  ! a negative initial value is a convention that indicates that groups (i.e. sub-communicators, one per run) are off by default
  integer :: mygroup = -1

  ! create a copy of the original output file path, to which we may add a "run0001/", "run0002/", "run0003/" prefix later
  ! if NUMBER_OF_SIMULTANEOUS_RUNS > 1
  character(len=MAX_STRING_LEN) :: OUTPUT_FILES = OUTPUT_FILES_BASE

  ! if doing simultaneous runs for the same mesh and model, see who should read the mesh and the model and broadcast it to others
  ! we put a default value here
  logical :: I_should_read_the_database = .true.

end module constants

!
!========================================================================
!

module shared_input_parameters

! holds input parameters given in DATA/Par_file

  use constants, only: MAX_STRING_LEN

  implicit none

  !--------------------------------------------------------------
  ! variables for preparation for compuation (simulation type, mesher, et al)
  !--------------------------------------------------------------
  ! name assigned to the running example
  character(len=MAX_STRING_LEN) :: title

  ! simulation type
  integer :: SIMULATION_TYPE

  ! NOISE_TOMOGRAPHY = 0 - turn noise tomography subroutines off; setting
  ! NOISE_TOMOGRAPHY equal to 0, in other words, results in an earthquake
  ! simulation rather than a noise simulation
  !
  ! NOISE_TOMOGRAPHY = 1 - compute "generating" wavefield and store the result;
  ! this stored wavefield is then used to compute the "ensemble forward"
  ! wavefield in the next noise simulation
  !
  ! NOISE_TOMOGRAPHY = 2 - compute "ensemble forward" wavefield and store the
  ! result; if an adjoint simulation is planned, users need to store
  ! seismograms (actually, cross-correlograms) for later processing
  !
  ! NOISE_TOMOGRAPHY = 3 - carry out adjoint simulation; users need to supply
  ! adjoint sources constructed from cross-correlograms computed during the
  ! "ensemble forward" step
  !
  ! For an explanation of terms and concepts in noise tomography, see "Tromp et
  ! al., 2011, Noise Cross-Correlation Sensitivity Kernels, Geophysical Journal
  ! International"
  integer :: NOISE_TOMOGRAPHY

  ! save forward arrays at the end of the simulation
  logical :: SAVE_FORWARD

  ! variables used for partitioning
  integer :: NPROC, partitioning_method

  ! number of control nodes
  integer :: ngnod

  ! number of time steps
  integer :: NSTEP

  ! time step size
  double precision :: DT

  ! value of time_stepping_scheme to decide which time scheme will be used
  ! 1 = Newmark (2nd order),
  ! 2 = LDDRK4-6 (4th-order 6-stage low storage Runge-Kutta),
  ! 3 = classical 4th-order 4-stage Runge-Kutta
  integer :: time_stepping_scheme

  ! simulation
  logical :: AXISYM

  logical :: P_SV

  ! computational platform type
  logical :: GPU_MODE

  ! creates/reads a binary database that allows to skip all time consuming setup steps in initialization
  ! 0 = does not read/create database
  ! 1 = creates database
  ! 2 = reads database
  integer :: setup_with_binary_database

  ! mesh files when using external mesh
  character(len=MAX_STRING_LEN) :: MODEL, SAVE_MODEL

  !#-----------------------------------------------------------------------------
  !#
  !# attenuation
  !#
  !#-----------------------------------------------------------------------------
  ! variables used for attenuation
  logical :: ATTENUATION_VISCOELASTIC
  logical :: ATTENUATION_PORO_FLUID_PART
  logical :: ATTENUATION_VISCOACOUSTIC
  double precision :: Q0_poroelastic,freq0_poroelastic

  integer :: N_SLS
  double precision :: ATTENUATION_f0_REFERENCE
  logical :: READ_VELOCITIES_AT_f0
  logical :: USE_SOLVOPT

  ! undo attenuation
  logical :: UNDO_ATTENUATION_AND_OR_PML
  ! variables used for iteration
  integer :: NT_DUMP_ATTENUATION

  !#-----------------------------------------------------------------------------
  !#
  !# sources
  !#
  !#-----------------------------------------------------------------------------
  ! variables used for source-receiver geometry
  integer :: NSOURCES
  logical :: force_normal_to_surface

  ! variables used for plane wave incidence
  logical :: initialfield
  logical :: add_Bielak_conditions_bottom,add_Bielak_conditions_right,add_Bielak_conditions_top,add_Bielak_conditions_left

  ! acoustic forcing of an acoustic medium at a rigid interface
  logical :: ACOUSTIC_FORCING

  !#-----------------------------------------------------------------------------
  !#
  !# receivers
  !#
  !#-----------------------------------------------------------------------------
  integer :: seismotype
  ! subsampling
  integer :: subsamp_seismos

  ! for better accuracy of pressure output (uses 2nd time-derivatives of the initial source time function)
  logical :: USE_TRICK_FOR_BETTER_PRESSURE

  integer :: NSTEP_BETWEEN_OUTPUT_SEISMOS

  ! Integrated energy field output
  logical :: COMPUTE_INTEGRATED_ENERGY_FIELD

  ! use this t0 as earliest starting time rather than the automatically calculated one
  ! (must be positive and bigger than the automatically one to be effective;
  !  simulation will start at t = - t0)
  double precision :: USER_T0

  ! seismogram format
  logical :: save_ASCII_seismograms
  logical :: save_binary_seismograms_single,save_binary_seismograms_double
  ! output seismograms in Seismic Unix format (adjoint traces will be read in the same format)
  logical :: SU_FORMAT

  logical :: use_existing_STATIONS

  integer :: nreceiversets

  double precision :: anglerec
  logical :: rec_normal_to_surface

  ! receiver sets
  integer, dimension(:),allocatable :: nrec_line
  double precision, dimension(:),allocatable :: xdeb,zdeb,xfin,zfin
  logical, dimension(:),allocatable :: record_at_surface_same_vertical

  !#-----------------------------------------------------------------------------
  !#
  !# adjoint kernel outputs
  !#
  !#-----------------------------------------------------------------------------
  ! kernel output in case of adjoint simulation
  logical :: save_ASCII_kernels

  integer :: NSTEP_BETWEEN_COMPUTE_KERNELS

  logical :: NO_BACKWARD_RECONSTRUCTION

  !#-----------------------------------------------------------------------------
  !#
  !# boundary conditions
  !#
  !#-----------------------------------------------------------------------------

  ! PML
  logical :: PML_BOUNDARY_CONDITIONS
  integer :: NELEM_PML_THICKNESS
  logical :: ROTATE_PML_ACTIVATE
  double precision :: ROTATE_PML_ANGLE
  double precision :: K_MIN_PML
  double precision :: K_MAX_PML
  double precision :: damping_change_factor_acoustic
  double precision :: damping_change_factor_elastic
  logical :: PML_PARAMETER_ADJUSTMENT

  ! Stacey
  logical :: STACEY_ABSORBING_CONDITIONS

  ! for horizontal periodic conditions: detect common points between left and right edges
  logical :: ADD_PERIODIC_CONDITIONS
  ! horizontal periodicity distance for periodic conditions
  double precision :: PERIODIC_HORIZ_DIST

  !#-----------------------------------------------------------------------------
  !#
  !# velocity and density models
  !#
  !#-----------------------------------------------------------------------------
  ! to store density and velocity model
  ! (actual material table will be read in in src/meshfem2D/read_material_table.f90)
  integer :: nbmodels

  ! input file name of TOMOGRAPHY
  character(len=MAX_STRING_LEN) :: TOMOGRAPHY_FILE

  logical :: read_external_mesh

  !#-----------------------------------------------------------------------------
  !#
  !# PARAMETERS FOR EXTERNAL MESHING
  !#
  !#-----------------------------------------------------------------------------
  character(len=MAX_STRING_LEN) :: mesh_file, nodes_coords_file, materials_file

  character(len=MAX_STRING_LEN) :: free_surface_file
  character(len=MAX_STRING_LEN) :: absorbing_surface_file
  character(len=MAX_STRING_LEN) :: acoustic_forcing_surface_file
  character(len=MAX_STRING_LEN) :: axial_elements_file
  character(len=MAX_STRING_LEN) :: absorbing_cpml_file
  character(len=MAX_STRING_LEN) :: tangential_detection_curve_file

  !#-----------------------------------------------------------------------------
  !#
  !# PARAMETERS FOR INTERNAL MESHING
  !#
  !#-----------------------------------------------------------------------------
  ! input parameter for in-house mesher
  character(len=MAX_STRING_LEN) :: interfacesfile

  double precision :: xmin_param,xmax_param
  integer :: nx_param

  ! variables used for absorbing boundary condition
  logical :: absorbbottom,absorbright,absorbtop,absorbleft

  ! number of regions
  ! (see reading in of regions table in read_regions.f90 file)
  integer :: nbregions

  !#-----------------------------------------------------------------------------
  !#
  !# display parameters
  !#
  !#-----------------------------------------------------------------------------
  ! general information during the computation and for information of the stability behavior during the simulation
  integer :: NSTEP_BETWEEN_OUTPUT_INFO

  ! for later check of the grid
  logical :: output_grid_Gnuplot,output_grid_ASCII

  ! for plotting the curve of energy
  logical :: OUTPUT_ENERGY
  integer :: NTSTEP_BETWEEN_OUTPUT_ENERGY

  !#-----------------------------------------------------------------------------
  !#
  !# movies/images/snaphots
  !#
  !#-----------------------------------------------------------------------------
  ! time step interval for image output
  integer :: NSTEP_BETWEEN_OUTPUT_IMAGES

  ! threshold value
  double precision :: cutsnaps

  ! JPEG image
  logical :: output_color_image
  integer :: imagetype_JPEG
  ! factor to subsample color images output by the code (useful for very large models)
  double precision :: factor_subsample_image
  ! by default the code normalizes each image independently to its maximum; use this option to use the global maximum below instead
  logical :: USE_CONSTANT_MAX_AMPLITUDE
  ! constant maximum amplitude to use for all color images if the USE_CONSTANT_MAX_AMPLITUDE option is true
  double precision :: CONSTANT_MAX_AMPLITUDE_TO_USE
  ! nonlinear display to enhance small amplitudes in color images
  double precision :: POWER_DISPLAY_COLOR
  logical :: DRAW_SOURCES_AND_RECEIVERS
  ! display acoustic layers as constant blue, because they likely correspond to water in the case of ocean acoustics
  ! or in the case of offshore oil industry experiments.
  ! (if off, display them as greyscale, as for elastic or poroelastic elements)
  logical :: DRAW_WATER_IN_BLUE
  ! use snapshot number in the file name of JPG color snapshots instead of the time step
  logical :: USE_SNAPSHOT_NUMBER_IN_FILENAME

  ! Postscript image
  logical :: output_postscript_snapshot
  integer :: imagetype_postscript
  logical :: meshvect,modelvect,boundvect,interpol
  ! number of interpolation points
  integer :: pointsdisp
  ! subsampling
  integer :: subsamp_postscript
  double precision :: sizemax_arrows
  ! US letter paper or European A4
  logical :: US_LETTER

  ! Wave field dumps
  logical :: output_wavefield_dumps
  integer :: imagetype_wavefield_dumps
  logical :: use_binary_for_wavefield_dumps

  ! NUMBER_OF_SIMULTANEOUS_RUNS
  integer :: NUMBER_OF_SIMULTANEOUS_RUNS
  logical :: BROADCAST_SAME_MESH_AND_MODEL

end module shared_input_parameters

!
!========================================================================
!

module shared_parameters

! for now, this just holds all input parameters given in DATA/Par_file
! in future, we might want both mesher and solver to read in the Par_file, similar to how 3D versions are handled

  use shared_input_parameters

  implicit none

  ! for MPI and partitioning
  integer :: myrank

  ! for Bielak condition
  logical :: add_Bielak_conditions

  ! for PML or Stacey boundary condition
  logical :: any_abs

  ! for interpolated snapshot
  logical :: plot_lowerleft_corner_only

  ! material file for changing the model parameter for inner mesh or updating the
  ! the material for an existed mesh
  ! (obsolete in Par_file now...)
  !logical :: assign_external_model, READ_EXTERNAL_SEP_FILE

  ! to store density and velocity model
  integer, dimension(:),allocatable :: num_material
  integer, dimension(:),allocatable :: icodemat

  double precision, dimension(:),allocatable :: rho_s_read
  double precision, dimension(:),allocatable :: rho_f_read

  ! acoustic/elastic/anisotropic
  double precision, dimension(:),allocatable :: cp,cs, &
    aniso3,aniso4,aniso5,aniso6,aniso7,aniso8,aniso9,aniso10,aniso11,aniso12,comp_g,QKappa,Qmu

  ! poroelastic
  ! note: adds ending _read to indicate these are readin values and to distinguish from solver arrays
  !       one could check if the solver arrays could be omitted and replaced with this ones in future...
  double precision, dimension(:),allocatable :: phi_read,tortuosity_read,permxx_read,permxz_read, &
       permzz_read,kappa_s_read,kappa_f_read,kappa_fr_read,eta_f_read,mu_fr_read

  ! mesh setup
  ! total number of elements
  integer :: nelmnts

  ! interface file data
  integer :: nx,nz
  integer :: nxread,nzread

  ! from interfaces file
  integer :: max_npoints_interface,number_of_interfaces

  ! vertical layers
  integer :: number_of_layers
  integer, dimension(:), allocatable :: nz_layer

end module shared_parameters


