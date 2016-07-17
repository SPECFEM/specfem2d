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


module constants

  include "constants.h"

end module constants

!
!========================================================================
!

module shared_input_parameters

! holds input parameters given in DATA/Par_file

  use constants,only: MAX_STRING_LEN

  implicit none

  !--------------------------------------------------------------
  ! variables for preparation for compuation (simulation type, mesher, et al)
  !--------------------------------------------------------------
  ! name assigned to the running example
  character(len=MAX_STRING_LEN) :: title

  ! simulation type
  integer :: SIMULATION_TYPE,NOISE_TOMOGRAPHY
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
  ! # 1 = Newmark (2nd order), 2 = LDDRK4-6 (4th-order 6-stage low storage Runge-Kutta)
  ! 3 = classical 4th-order 4-stage Runge-Kutta
  integer :: time_stepping_scheme

  ! simulation
  logical :: AXISYM

  logical :: P_SV

  ! computational platform type
  logical :: GPU_MODE

  ! mesh files when using external mesh
  character(len=MAX_STRING_LEN) :: MODEL, SAVE_MODEL

  !#-----------------------------------------------------------------------------
  !#
  !# attenuation
  !#
  !#-----------------------------------------------------------------------------
  ! variables used for attenuation
  logical :: ATTENUATION_VISCOELASTIC_SOLID
  logical :: ATTENUATION_PORO_FLUID_PART
  double precision :: Q0,freq0

  integer :: N_SLS
  double precision :: f0_attenuation
  logical :: READ_VELOCITIES_AT_f0

  ! undo attenuation
  logical :: UNDO_ATTENUATION
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
  integer, dimension(:),allocatable :: nrec
  double precision, dimension(:),allocatable :: xdeb,zdeb,xfin,zfin
  logical, dimension(:),allocatable :: record_at_surface_same_vertical

  !#-----------------------------------------------------------------------------
  !#
  !# adjoint kernel outputs
  !#
  !#-----------------------------------------------------------------------------
  ! kernel output in case of adjoint simulation
  logical :: save_ASCII_kernels


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
  ! (actual material table will be read in in src/meshfem2D/read_materials.f90)
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
  character(len=MAX_STRING_LEN) :: CPML_element_file
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
  ! general information during the computation
  integer :: NSTEP_BETWEEN_OUTPUT_INFO

  ! for later check of the grid
  logical :: output_grid_Gnuplot,output_grid_ASCII

  ! for plotting the curve of energy
  logical :: output_energy

  !#-----------------------------------------------------------------------------
  !#
  !# movies/images/snaphots
  !#
  !#-----------------------------------------------------------------------------

  integer :: NSTEP_BETWEEN_OUTPUT_IMAGES

  integer :: NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS

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
  ! non linear display to enhance small amplitudes in color images
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
  integer :: pointsdisp,subsamp_postscript
  double precision :: sizemax_arrows
  ! US letter paper or European A4
  logical :: US_LETTER

  ! Wave field dumps
  logical :: output_wavefield_dumps
  integer :: imagetype_wavefield_dumps
  logical :: use_binary_for_wavefield_dumps

end module shared_input_parameters

!
!========================================================================
!

module shared_parameters

! for now, this just holds all input parameters given in DATA/Par_file
! in future, we might want both mesher and solver to read in the Par_file, similar to how 3D versions are handled

  use shared_input_parameters

  implicit none

end module shared_parameters


