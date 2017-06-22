!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
! Copyright CNRS, Inria and University of Pau, France,
! and Princeton University / California Institute of Technology, USA.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
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

module parameter_file

  ! note: we use this module definition only to be able to allocate
  !          arrays for receiverlines and materials in this subroutine rather than in the main
  !          routine in meshfem2D.F90

  ! note 2: the filename ending is .F90 to have pre-compilation with pragmas
  !            (like #ifndef USE_MPI) working properly

  implicit none
  character(len=100) :: interfacesfile,title

  integer :: SIMULATION_TYPE, NOISE_TOMOGRAPHY
  logical :: SAVE_FORWARD,read_external_mesh

  character(len=256) :: mesh_file, nodes_coords_file, materials_file, &
                        free_surface_file, absorbing_surface_file, &
                        CPML_element_file
  character(len=256)  :: tangential_detection_curve_file

  ! variables used for partitioning
  integer :: nproc,partitioning_method

  double precision :: xmin,xmax
  integer :: nx,ngnod

  logical :: initialfield,add_Bielak_conditions,assign_external_model, &
            READ_EXTERNAL_SEP_FILE,ATTENUATION_VISCOELASTIC_SOLID,ATTENUATION_PORO_FLUID_PART, &
            save_ASCII_seismograms,save_binary_seismograms_single,save_binary_seismograms_double,DRAW_SOURCES_AND_RECEIVERS

  double precision :: Q0,freq0

  logical :: P_SV
  logical :: any_abs,absbottom,absright,abstop,absleft

  integer :: nt
  double precision :: deltat

  integer :: NSOURCES
  logical :: force_normal_to_surface

  ! variables used for attenuation
  integer  :: N_SLS
  double precision  :: f0_attenuation

  integer :: seismotype
  logical :: seismo_p
  logical :: generate_STATIONS

  integer :: nreceiversets
  double precision :: anglerec
  logical :: rec_normal_to_surface

  integer, dimension(:), pointer :: nrec
  double precision, dimension(:), pointer :: xdeb,zdeb,xfin,zfin
  logical, dimension(:), pointer :: enreg_surf_same_vertical

  integer :: NSTEP_BETWEEN_OUTPUT_INFO,NSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP_BETWEEN_OUTPUT_IMAGES,NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS, &
             subsamp_seismos,imagetype_JPEG,imagetype_wavefield_dumps,NELEM_PML_THICKNESS
  logical :: output_postscript_snapshot,output_color_image,PML_BOUNDARY_CONDITIONS
  integer :: imagetype_postscript
  double precision :: cutsnaps
  logical :: meshvect,modelvect,boundvect,interpol
  integer :: pointsdisp,subsamp_postscript
  double precision :: sizemax_arrows
  logical :: output_grid_Gnuplot,output_grid_ASCII,output_energy,output_wavefield_dumps,use_binary_for_wavefield_dumps
  logical :: plot_lowerleft_corner_only

  ! to store density and velocity model
  integer :: nb_materials
  integer, dimension(:),pointer :: icodemat
  double precision, dimension(:),pointer :: rho_s,cp,cs, &
    aniso3,aniso4,aniso5,aniso6,aniso7,aniso8,QKappa,Qmu
  double precision, dimension(:),pointer :: rho_f,phi,tortuosity,permxx,permxz, &
       permzz,kappa_s,kappa_f,kappa_fr,eta_f,mu_fr

! factor to subsample color images output by the code (useful for very large models)
  integer :: factor_subsample_image

! use snapshot number in the file name of JPG color snapshots instead of the time step
  logical :: USE_SNAPSHOT_NUMBER_IN_FILENAME

! display acoustic layers as constant blue, because they likely correspond to water in the case of ocean acoustics
! or in the case of offshore oil industry experiments.
! (if off, display them as greyscale, as for elastic or poroelastic elements)
  logical :: DRAW_WATER_IN_BLUE

! US letter paper or European A4
  logical :: US_LETTER

! non linear display to enhance small amplitudes in color images
  double precision :: POWER_DISPLAY_COLOR

! perform inverse Cuthill-McKee (1969) permutation for mesh numbering
  logical :: PERFORM_CUTHILL_MCKEE

! output seismograms in Seismic Unix format (adjoint traces will be read in the same format)
  logical :: SU_FORMAT

! use this t0 as earliest starting time rather than the automatically calculated one
! (must be positive and bigger than the automatically one to be effective;
!  simulation will start at t = - t0)
  double precision :: USER_T0

! value of time_stepping_scheme to decide which time scheme will be used
! # 1 = Newmark (2nd order), 2 = LDDRK4-6 (4th-order 6-stage low storage Runge-Kutta)
! 3 = classical 4th-order 4-stage Runge-Kutta
  integer :: time_stepping_scheme

!! DK DK for horizontal periodic conditions: detect common points between left and right edges
  logical :: ADD_PERIODIC_CONDITIONS

!! DK DK horizontal periodicity distance for periodic conditions
  double precision :: PERIODIC_horiz_dist

!! DK DK grid point detection tolerance for periodic conditions
  double precision :: PERIODIC_DETECT_TOL

contains

  subroutine read_parameter_file()

! reads in DATA/Par_file

  implicit none
  include "constants.h"

  ! local parameters
  integer :: ios,ireceiverlines
  integer,external :: err_occurred

  ! read file names and path for output
  call read_value_string_p(title, 'solver.title')
  if (err_occurred() /= 0) stop 'error reading parameter 1 in Par_file'

  write(*,*) 'Title of the simulation'
  write(*,*) title
  print *

  ! read type of simulation
  call read_value_integer_p(SIMULATION_TYPE, 'solver.SIMULATION_TYPE')
  if (err_occurred() /= 0) stop 'error reading parameter 2 in Par_file'

  call read_value_integer_p(NOISE_TOMOGRAPHY, 'solver.NOISE_TOMOGRAPHY')
  if (err_occurred() /= 0) stop 'error reading parameter NOISE_TOMOGRAPHY in Par_file'

  call read_value_logical_p(SAVE_FORWARD, 'solver.SAVE_FORWARD')
  if (err_occurred() /= 0) stop 'error reading parameter 3 in Par_file'

  ! read info about partitioning
  call read_value_integer_p(nproc, 'solver.nproc')
  if (err_occurred() /= 0) stop 'error reading parameter 4 in Par_file'

  call read_value_integer_p(partitioning_method, 'mesher.partitioning_method')
  if (err_occurred() /= 0) stop 'error reading parameter 5a in Par_file'

  call read_value_logical_p(PERFORM_CUTHILL_MCKEE, 'mesher.PERFORM_CUTHILL_MCKEE')
  if (err_occurred() /= 0) stop 'error reading parameter 5b in Par_file'

  call read_value_integer_p(ngnod, 'mesher.ngnod')
  if (err_occurred() /= 0) stop 'error reading parameter 6 in Par_file'

  call read_value_logical_p(initialfield, 'solver.initialfield')
  if (err_occurred() /= 0) stop 'error reading parameter 7 in Par_file'

  call read_value_logical_p(add_Bielak_conditions, 'solver.add_Bielak_conditions')
  if (err_occurred() /= 0) stop 'error reading parameter 8 in Par_file'

  call read_value_logical_p(assign_external_model, 'mesher.assign_external_model')
  if (err_occurred() /= 0) stop 'error reading parameter 9 in Par_file'

  call read_value_logical_p(READ_EXTERNAL_SEP_FILE, 'mesher.READ_EXTERNAL_SEP_FILE')
  if (err_occurred() /= 0) stop 'error reading parameter 10 in Par_file'

  call read_value_logical_p(ATTENUATION_VISCOELASTIC_SOLID, 'solver.ATTENUATION_VISCOELASTIC_SOLID')
  if (err_occurred() /= 0) stop 'error reading parameter 11 in Par_file'

  ! read viscous attenuation parameters (poroelastic media)
  call read_value_logical_p(ATTENUATION_PORO_FLUID_PART, 'solver.ATTENUATION_PORO_FLUID_PART')
  if (err_occurred() /= 0) stop 'error reading parameter 12a in Par_file'

  call read_value_double_precision_p(Q0, 'solver.Q0')
  if (err_occurred() /= 0) stop 'error reading parameter 13 in Par_file'

  call read_value_double_precision_p(freq0, 'solver.freq0')
  if (err_occurred() /= 0) stop 'error reading parameter 14 in Par_file'

  ! determine if body or surface (membrane) waves calculation
  call read_value_logical_p(P_SV, 'solver.P_SV')
  if (err_occurred() /= 0) stop 'error reading parameter 15 in Par_file'

  ! read time step parameters
  call read_value_integer_p(nt, 'solver.nt')
  if (err_occurred() /= 0) stop 'error reading parameter 16 in Par_file'

  call read_value_double_precision_p(deltat, 'solver.deltat')
  if (err_occurred() /= 0) stop 'error reading parameter 17a in Par_file'

  call read_value_double_precision_p(USER_T0, 'solver.USER_T0')
  if (err_occurred() /= 0) stop 'error reading parameter 17b in Par_file'

  call read_value_integer_p(time_stepping_scheme, 'solver.time_stepping_scheme')
  if (err_occurred() /= 0) stop 'error reading parameter 17c in Par_file'

  ! read source infos
  call read_value_integer_p(NSOURCES, 'solver.NSOURCES')
  if (err_occurred() /= 0) stop 'error reading parameter 18 in Par_file'

  call read_value_logical_p(force_normal_to_surface, 'solver.force_normal_to_surface')
  if (err_occurred() /= 0) stop 'error reading parameter 19 in Par_file'

  ! read constants for attenuation
  call read_value_integer_p(N_SLS, 'solver.N_SLS')
  if (err_occurred() /= 0) stop 'error reading parameter 20 in Par_file'
  if (N_SLS < 1) stop 'must have N_SLS >= 1 even if attenuation if off because it is used to assign some arrays'

  call read_value_double_precision_p(f0_attenuation, 'solver.f0_attenuation')
  if (err_occurred() /= 0) stop 'error reading parameter 21 in Par_file'

  ! read receiver line parameters
  call read_value_integer_p(seismotype, 'solver.seismotype')
  if (err_occurred() /= 0) stop 'error reading parameter 22 in Par_file'

  !RC: Read whether pressure is also recorded at each receiver. Note that
  !this will fail for "old" parameter files that don't have seismo_p.
  call read_value_logical_p(seismo_p, 'solver.seismo_p')
  if (err_occurred() /= 0) stop 'error reading parameter 22b in Par_file'

  call read_value_integer_p(NSTEP_BETWEEN_OUTPUT_SEISMOS, 'solver.NSTEP_BETWEEN_OUTPUT_SEISMOS')
  if (err_occurred() /= 0) stop 'error reading parameter 33b in Par_file'

  call read_value_logical_p(save_ASCII_seismograms, 'solver.save_ASCII_seismograms')
  if (err_occurred() /= 0) stop 'error reading parameter 12b in Par_file'

  call read_value_logical_p(save_binary_seismograms_single, 'solver.save_binary_seismograms_single')
  if (err_occurred() /= 0) stop 'error reading parameter 12c in Par_file'

  call read_value_logical_p(save_binary_seismograms_double, 'solver.save_binary_seismograms_double')
  if (err_occurred() /= 0) stop 'error reading parameter 12cc in Par_file'

  call read_value_logical_p(SU_FORMAT, 'solver.SU_FORMAT')
  if (err_occurred() /= 0) stop 'error reading parameter 26b in Par_file'

  call read_value_integer_p(subsamp_seismos, 'solver.subsamp_seismos')
  if (err_occurred() /= 0) stop 'error reading parameter 33e in Par_file'
  if (subsamp_seismos < 1) stop 'error: subsamp_seismos must be >= 1'

  call read_value_logical_p(generate_STATIONS, 'solver.generate_STATIONS')
  if (err_occurred() /= 0) stop 'error reading parameter 23 in Par_file'

  call read_value_integer_p(nreceiversets, 'solver.nreceiversets')
  if (err_occurred() /= 0) stop 'error reading parameter 24 in Par_file'

  call read_value_double_precision_p(anglerec, 'solver.anglerec')
  if (err_occurred() /= 0) stop 'error reading parameter 25 in Par_file'

  call read_value_logical_p(rec_normal_to_surface, 'solver.rec_normal_to_surface')
  if (err_occurred() /= 0) stop 'error reading parameter 26a in Par_file'

  if (nreceiversets < 1) stop 'number of receiver lines must be greater than 1'

  ! allocate receiver line arrays
  allocate(nrec(nreceiversets))
  allocate(xdeb(nreceiversets))
  allocate(zdeb(nreceiversets))
  allocate(xfin(nreceiversets))
  allocate(zfin(nreceiversets))
  allocate(enreg_surf_same_vertical(nreceiversets),stat=ios)
  if ( ios /= 0 ) stop 'error allocating receiver lines'

  ! loop on all the receiver lines
  do ireceiverlines = 1,nreceiversets
    call read_value_integer_next_p(nrec(ireceiverlines),'solver.nrec')
    if (err_occurred() /= 0) stop 'error reading parameter 27 in Par_file'

    call read_value_double_prec_next_p(xdeb(ireceiverlines),'solver.xdeb')
    if (err_occurred() /= 0) stop 'error reading parameter 28 in Par_file'

    call read_value_double_prec_next_p(zdeb(ireceiverlines),'solver.zdeb')
    if (err_occurred() /= 0) stop 'error reading parameter 29 in Par_file'

    call read_value_double_prec_next_p(xfin(ireceiverlines),'solver.xfin')
    if (err_occurred() /= 0) stop 'error reading parameter 30 in Par_file'

    call read_value_double_prec_next_p(zfin(ireceiverlines),'solver.zfin')
    if (err_occurred() /= 0) stop 'error reading parameter 31 in Par_file'

    call read_value_logical_next_p(enreg_surf_same_vertical(ireceiverlines),'solver.enreg_surf_same_vertical')
    if (err_occurred() /= 0) stop 'error reading parameter 32 in Par_file'

    if (read_external_mesh .and. enreg_surf_same_vertical(ireceiverlines)) then
      stop 'Cannot use enreg_surf_same_vertical with external meshes!'
    endif
  enddo

  ! read display parameters
  call read_value_integer_p(NSTEP_BETWEEN_OUTPUT_INFO, 'solver.NSTEP_BETWEEN_OUTPUT_INFO')
  if (err_occurred() /= 0) stop 'error reading parameter 33a in Par_file'

  call read_value_integer_p(NSTEP_BETWEEN_OUTPUT_IMAGES, 'solver.NSTEP_BETWEEN_OUTPUT_IMAGES')
  if (err_occurred() /= 0) stop 'error reading parameter 33c in Par_file'

  call read_value_double_precision_p(cutsnaps, 'solver.cutsnaps')
  if (err_occurred() /= 0) stop 'error reading parameter 37 in Par_file'

  call read_value_logical_p(output_color_image, 'solver.output_color_image')
  if (err_occurred() /= 0) stop 'error reading parameter 35 in Par_file'

  call read_value_integer_p(imagetype_JPEG, 'solver.imagetype_JPEG')
  if (err_occurred() /= 0) stop 'error reading parameter 33f in Par_file'

  call read_value_integer_p(factor_subsample_image, 'solver.factor_subsample_image')
  if (err_occurred() /= 0) stop 'error reading parameter 43b in Par_file'

  call read_value_double_precision_p(POWER_DISPLAY_COLOR, 'solver.POWER_DISPLAY_COLOR')
  if (err_occurred() /= 0) stop 'error reading parameter 43c in Par_file'

  call read_value_logical_p(DRAW_SOURCES_AND_RECEIVERS, 'solver.DRAW_SOURCES_AND_RECEIVERS')
  if (err_occurred() /= 0) stop 'error reading parameter 12d in Par_file'

  call read_value_logical_p(DRAW_WATER_IN_BLUE, 'solver.DRAW_WATER_IN_BLUE')
  if (err_occurred() /= 0) stop 'error reading parameter 43d in Par_file'

  call read_value_logical_p(USE_SNAPSHOT_NUMBER_IN_FILENAME, 'solver.USE_SNAPSHOT_NUMBER_IN_FILENAME')
  if (err_occurred() /= 0) stop 'error reading parameter 44c in Par_file'

  call read_value_logical_p(output_postscript_snapshot, 'solver.output_postscript_snapshot')
  if (err_occurred() /= 0) stop 'error reading parameter 34 in Par_file'

  call read_value_integer_p(imagetype_postscript, 'solver.imagetype_postscript')
  if (err_occurred() /= 0) stop 'error reading parameter 36 in Par_file'

  call read_value_logical_p(meshvect, 'solver.meshvect')
  if (err_occurred() /= 0) stop 'error reading parameter 38 in Par_file'

  call read_value_logical_p(modelvect, 'solver.modelvect')
  if (err_occurred() /= 0) stop 'error reading parameter 39 in Par_file'

  call read_value_logical_p(boundvect, 'solver.boundvect')
  if (err_occurred() /= 0) stop 'error reading parameter 40 in Par_file'

  call read_value_logical_p(interpol, 'solver.interpol')
  if (err_occurred() /= 0) stop 'error reading parameter 41 in Par_file'

  call read_value_integer_p(pointsdisp, 'solver.pointsdisp')
  if (err_occurred() /= 0) stop 'error reading parameter 42 in Par_file'

  call read_value_integer_p(subsamp_postscript, 'solver.subsamp_postscript')
  if (err_occurred() /= 0) stop 'error reading parameter 43a in Par_file'

  call read_value_double_precision_p(sizemax_arrows, 'solver.sizemax_arrows')
  if (err_occurred() /= 0) stop 'error reading parameter 44a in Par_file'

  call read_value_logical_p(US_LETTER, 'solver.US_LETTER')
  if (err_occurred() /= 0) stop 'error reading parameter 44b in Par_file'

  call read_value_integer_p(NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS, 'solver.NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS')
  if (err_occurred() /= 0) stop 'error reading parameter 33d in Par_file'

  call read_value_logical_p(output_wavefield_dumps, 'solver.output_wavefield_dumps')
  if (err_occurred() /= 0) stop 'error reading parameter 48 in Par_file'

  call read_value_integer_p(imagetype_wavefield_dumps, 'solver.imagetype_wavefield_dumps')
  if (err_occurred() /= 0) stop 'error reading parameter 33g in Par_file'

  call read_value_logical_p(use_binary_for_wavefield_dumps, 'solver.use_binary_for_wavefield_dumps')
  if (err_occurred() /= 0) stop 'error reading parameter 48 in Par_file'

  call read_value_logical_p(output_grid_Gnuplot, 'solver.output_grid_Gnuplot')
  if (err_occurred() /= 0) stop 'error reading parameter 45 in Par_file'

  call read_value_logical_p(output_grid_ASCII, 'solver.output_grid_ASCII')
  if (err_occurred() /= 0) stop 'error reading parameter 46 in Par_file'

  call read_value_logical_p(output_energy, 'solver.output_energy')
  if (err_occurred() /= 0) stop 'error reading parameter 47 in Par_file'

  ! read the different material materials
  call read_value_integer_p(nb_materials, 'mesher.nbmodels')
  if (err_occurred() /= 0) stop 'error reading parameter 49 in Par_file'

  if (nb_materials <= 0) stop 'Non-positive number of materials not allowed!'

  allocate(icodemat(nb_materials))
  allocate(cp(nb_materials))
  allocate(cs(nb_materials))
  allocate(aniso3(nb_materials))
  allocate(aniso4(nb_materials))
  allocate(aniso5(nb_materials))
  allocate(aniso6(nb_materials))
  allocate(aniso7(nb_materials))
  allocate(aniso8(nb_materials))
  allocate(QKappa(nb_materials))
  allocate(Qmu(nb_materials))
  allocate(rho_s(nb_materials))
  allocate(rho_f(nb_materials))
  allocate(phi(nb_materials))
  allocate(tortuosity(nb_materials))
  allocate(permxx(nb_materials))
  allocate(permxz(nb_materials))
  allocate(permzz(nb_materials))
  allocate(kappa_s(nb_materials))
  allocate(kappa_f(nb_materials))
  allocate(kappa_fr(nb_materials))
  allocate(eta_f(nb_materials))
  allocate(mu_fr(nb_materials))

  call read_materials(nb_materials,icodemat,cp,cs, &
                      aniso3,aniso4,aniso5,aniso6,aniso7,aniso8, &
                      QKappa,Qmu,rho_s,rho_f,phi,tortuosity, &
                      permxx,permxz,permzz,kappa_s,kappa_f,kappa_fr, &
                      eta_f,mu_fr)

  ! boolean defining whether internal or external mesh
  call read_value_logical_p(read_external_mesh, 'mesher.read_external_mesh')
  if (err_occurred() /= 0) stop 'error reading parameter 50 in Par_file'

  call read_value_logical_p(PML_BOUNDARY_CONDITIONS, 'solver.PML_BOUNDARY_CONDITIONS')
  if (err_occurred() /= 0) stop 'error reading parameter 33za in Par_file'

  call read_value_integer_p(NELEM_PML_THICKNESS, 'solver.NELEM_PML_THICKNESS')
  if (err_occurred() /= 0) stop 'error reading parameter 33zb in Par_file'

  ! boolean defining whether to use any absorbing boundaries
  call read_value_logical_p(any_abs, 'solver.STACEY_ABSORBING_CONDITIONS')
  if (err_occurred() /= 0) stop 'error reading parameter 51a in Par_file'

  if (add_Bielak_conditions .and. .not. any_abs) &
    stop 'need STACEY_ABSORBING_CONDITIONS set to .true. in order to use add_Bielak_conditions'

  ! solve the conflict in value of PML_BOUNDARY_CONDITIONS and STACEY_ABSORBING_CONDITIONS
  if (PML_BOUNDARY_CONDITIONS) any_abs = .true.

  call read_value_logical_p(ADD_PERIODIC_CONDITIONS, 'solver.ADD_PERIODIC_CONDITIONS')
  if (err_occurred() /= 0) stop 'error reading parameter 51b in Par_file'

  call read_value_double_precision_p(PERIODIC_horiz_dist, 'solver.PERIODIC_horiz_dist')
  if (err_occurred() /= 0) stop 'error reading parameter 51c in Par_file'

  call read_value_double_precision_p(PERIODIC_DETECT_TOL, 'solver.PERIODIC_DETECT_TOL')
  if (err_occurred() /= 0) stop 'error reading parameter 51d in Par_file'

  !-----------------
  ! external mesh parameters

  if ( read_external_mesh ) then

  ! read info about external mesh
  call read_value_string_p(mesh_file, 'mesher.mesh_file')
  if (err_occurred() /= 0) stop 'error reading parameter 52 in Par_file'

  call read_value_string_p(nodes_coords_file, 'mesher.nodes_coords_file')
  if (err_occurred() /= 0) stop 'error reading parameter 53 in Par_file'

  call read_value_string_p(materials_file, 'mesher.materials_file')
  if (err_occurred() /= 0) stop 'error reading parameter 54 in Par_file'

  call read_value_string_p(free_surface_file, 'mesher.free_surface_file')
  if (err_occurred() /= 0) stop 'error reading parameter 55 in Par_file'

  call read_value_string_p(absorbing_surface_file, 'mesher.absorbing_surface_file')
  if (err_occurred() /= 0) stop 'error reading parameter 56 in Par_file'

  call read_value_string_p(CPML_element_file, 'mesher.CPML_element_file')
  if (err_occurred() /= 0) stop 'error reading parameter 56 in Par_file'

  call read_value_string_p(tangential_detection_curve_file, 'mesher.tangential_detection_curve_file')
  if (err_occurred() /= 0) stop 'error reading parameter 57 in Par_file'

  else

  !-----------------
  ! internal mesh parameters

  ! interfaces file
  call read_value_string_p(interfacesfile, 'mesher.interfacesfile')
  if (err_occurred() /= 0) stop 'error reading parameter 58 in Par_file'

  ! read grid parameters
  call read_value_double_precision_p(xmin, 'mesher.xmin')
  if (err_occurred() /= 0) stop 'error reading parameter 59 in Par_file'

  call read_value_double_precision_p(xmax, 'mesher.xmax')
  if (err_occurred() /= 0) stop 'error reading parameter 60 in Par_file'

  call read_value_integer_p(nx, 'mesher.nx')
  if (err_occurred() /= 0) stop 'error reading parameter 61 in Par_file'

  ! read absorbing boundary parameters
  call read_value_logical_p(absbottom, 'solver.absorbbottom')
  if (err_occurred() /= 0) stop 'error reading parameter 62 in Par_file'

  call read_value_logical_p(absright, 'solver.absorbright')
  if (err_occurred() /= 0) stop 'error reading parameter 63 in Par_file'

  call read_value_logical_p(abstop, 'solver.absorbtop')
  if (err_occurred() /= 0) stop 'error reading parameter 64 in Par_file'

  call read_value_logical_p(absleft, 'solver.absorbleft')
  if (err_occurred() /= 0) stop 'error reading parameter 65 in Par_file'

  ! note: if internal mesh, then regions will be read in by read_regions (from meshfem2D)

  endif

  ! checks input parameters
  call check_parameters()

  end subroutine read_parameter_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine check_parameters()

  implicit none

  ! checks partitioning
  if ( nproc <= 0 ) then
     print *, 'Number of processes (nproc) must be greater than or equal to one.'
     stop
  endif

#ifndef USE_MPI
  if ( nproc > 1 ) then
     print *, 'Number of processes (nproc) must be equal to one when not using MPI.'
     print *, 'Please recompile with -DUSE_MPI in order to enable use of MPI.'
     stop
  endif
#endif

  if (partitioning_method /= 1 .and. partitioning_method /= 3) then
     print *, 'Invalid partitioning method number.'
     print *, 'Partitioning method ',partitioning_method,' was requested, but is not available.'
     print *, 'Support for the METIS graph partitioner has been discontinued, please use SCOTCH (option 3) instead.'
     stop
  endif

  ! checks absorbing boundaries
  if (.not. any_abs ) then
     absbottom = .false.
     absright = .false.
     abstop = .false.
     absleft = .false.
  endif

  ! can use only one point to display lower-left corner only for interpolated snapshot
  if (pointsdisp < 3) then
     pointsdisp = 3
     plot_lowerleft_corner_only = .true.
  else
     plot_lowerleft_corner_only = .false.
  endif

  end subroutine check_parameters

end module parameter_file

