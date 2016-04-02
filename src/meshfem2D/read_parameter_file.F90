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

  subroutine read_parameter_file()

! reads in DATA/Par_file

  use parameter_file_par

  implicit none

  include "constants.h"

  ! local parameters
  integer :: ios,ireceiverlines
  integer,external :: err_occurred

  !--------------------------------------------------------------------
  !
  ! simulation input paramters
  !
  !--------------------------------------------------------------------

  ! read file names and path for output
  call read_value_string_p(title, 'solver.title')
  if (err_occurred() /= 0) stop 'error reading parameter title in Par_file'

  ! read type of simulation
  call read_value_integer_p(SIMULATION_TYPE, 'solver.SIMULATION_TYPE')
  if (err_occurred() /= 0) stop 'error reading parameter SIMULATION_TYPE in Par_file'

  call read_value_integer_p(NOISE_TOMOGRAPHY, 'solver.NOISE_TOMOGRAPHY')
  if (err_occurred() /= 0) stop 'error reading parameter NOISE_TOMOGRAPHY in Par_file'

  call read_value_logical_p(SAVE_FORWARD, 'solver.SAVE_FORWARD')
  if (err_occurred() /= 0) stop 'error reading parameter SAVE_FORWARD in Par_file'

  ! read info about partitioning
  call read_value_integer_p(NPROC, 'solver.NPROC')
  if (err_occurred() /= 0) stop 'error reading parameter NPROC in Par_file'

  call read_value_integer_p(partitioning_method, 'mesher.partitioning_method')
  if (err_occurred() /= 0) stop 'error reading parameter partitioning_method in Par_file'

  call read_value_integer_p(ngnod, 'mesher.ngnod')
  if (err_occurred() /= 0) stop 'error reading parameter ngnod in Par_file'

  ! read time step parameters
  call read_value_integer_p(NSTEP, 'solver.NSTEP')
  if (err_occurred() /= 0) stop 'error reading parameter NSTEP in Par_file'

  call read_value_double_precision_p(DT, 'solver.DT')
  if (err_occurred() /= 0) stop 'error reading parameter DT in Par_file'

  call read_value_integer_p(time_stepping_scheme, 'solver.time_stepping_scheme')
  if (err_occurred() /= 0) stop 'error reading parameter time_stepping_scheme in Par_file'

  ! axisymmetric (2.5D) or Cartesian planar (2D) simulation
  call read_value_logical_p(AXISYM, 'solver.AXISYM')
  if (err_occurred() /= 0) stop 'error reading parameter AXISYM in Par_file'

  ! determine if body or surface (membrane) waves calculation
  call read_value_logical_p(P_SV, 'solver.P_SV')
  if (err_occurred() /= 0) stop 'error reading parameter P_SV in Par_file'

  call read_value_logical_p(GPU_MODE, 'solver.GPU_MODE')
  if (err_occurred() /= 0) stop 'error reading parameter GPU_MODE in Par_file'

  call read_value_string_p(MODEL, 'mesher.MODEL')
  if (err_occurred() /= 0) stop 'error reading parameter MODEL in Par_file'

  call read_value_string_p(SAVE_MODEL, 'mesher.SAVE_MODEL')
  if (err_occurred() /= 0) stop 'error reading parameter SAVE_MODEL in Par_file'

  ! user output
  write(*,*) 'Title of the simulation'
  write(*,*) title
  print *
  if (AXISYM) write(*,*) 'Axisymmetric simulation'
  print *

  !--------------------------------------------------------------------
  !
  ! attenuation
  !
  !--------------------------------------------------------------------

  call read_value_logical_p(ATTENUATION_VISCOELASTIC_SOLID, 'solver.ATTENUATION_VISCOELASTIC_SOLID')
  if (err_occurred() /= 0) stop 'error reading parameter ATTENUATION_VISCOELASTIC_SOLID in Par_file'

  ! read viscous attenuation parameters (poroelastic media)
  call read_value_logical_p(ATTENUATION_PORO_FLUID_PART, 'solver.ATTENUATION_PORO_FLUID_PART')
  if (err_occurred() /= 0) stop 'error reading parameter ATTENUATION_PORO_FLUID_PART in Par_file'

  call read_value_double_precision_p(Q0, 'solver.Q0')
  if (err_occurred() /= 0) stop 'error reading parameter Q0 in Par_file'

  call read_value_double_precision_p(freq0, 'solver.freq0')
  if (err_occurred() /= 0) stop 'error reading parameter freq0 in Par_file'

  ! read constants for attenuation
  call read_value_integer_p(N_SLS, 'solver.N_SLS')
  if (err_occurred() /= 0) stop 'error reading parameter N_SLS in Par_file'

  call read_value_double_precision_p(f0_attenuation, 'solver.f0_attenuation')
  if (err_occurred() /= 0) stop 'error reading parameter f0_attenuation in Par_file'

  call read_value_logical_p(READ_VELOCITIES_AT_f0, 'solver.READ_VELOCITIES_AT_f0')
  if (err_occurred() /= 0) stop 'error reading parameter READ_VELOCITIES_AT_f0 in Par_file'

  call read_value_logical_p(UNDO_ATTENUATION, 'solver.UNDO_ATTENUATION')
  if (err_occurred() /= 0) stop 'error reading parameter UNDO_ATTENUATION in Par_file'

  call read_value_integer_p(NT_DUMP_ATTENUATION, 'solver.NT_DUMP_ATTENUATION')
  if (err_occurred() /= 0) stop 'error reading parameter NT_DUMP_ATTENUATION in Par_file'

  !--------------------------------------------------------------------
  !
  ! sources
  !
  !--------------------------------------------------------------------

  ! read source infos
  call read_value_integer_p(NSOURCES, 'solver.NSOURCES')
  if (err_occurred() /= 0) stop 'error reading parameter NSOURCES in Par_file'

  call read_value_logical_p(force_normal_to_surface, 'solver.force_normal_to_surface')
  if (err_occurred() /= 0) stop 'error reading parameter force_normal_to_surface in Par_file'

  call read_value_logical_p(initialfield, 'solver.initialfield')
  if (err_occurred() /= 0) stop 'error reading parameter initialfield in Par_file'
  if (initialfield .and. NPROC > 1) stop 'initialfield (plane waves) currently have a bug in parallel i.e. when &
     & NPROC > 1, see https://github.com/geodynamics/specfem2d/issues/550 , thus please run with a single processor'

  call read_value_logical_p(add_Bielak_conditions_bottom, 'solver.add_Bielak_conditions_bottom')
  if (err_occurred() /= 0) stop 'error reading parameter add_Bielak_conditions_bottom in Par_file'

  call read_value_logical_p(add_Bielak_conditions_right, 'solver.add_Bielak_conditions_right')
  if (err_occurred() /= 0) stop 'error reading parameter add_Bielak_conditions_right in Par_file'

  call read_value_logical_p(add_Bielak_conditions_top, 'solver.add_Bielak_conditions_top')
  if (err_occurred() /= 0) stop 'error reading parameter add_Bielak_conditions_top in Par_file'

  call read_value_logical_p(add_Bielak_conditions_left, 'solver.add_Bielak_conditions_left')
  if (err_occurred() /= 0) stop 'error reading parameter add_Bielak_conditions_left in Par_file'

  add_Bielak_conditions = add_Bielak_conditions_bottom .or. add_Bielak_conditions_right .or. &
                          add_Bielak_conditions_top .or. add_Bielak_conditions_left

  ! read acoustic forcing flag
  call read_value_logical_p(ACOUSTIC_FORCING, 'solver.ACOUSTIC_FORCING')
  if (err_occurred() /= 0) stop 'error reading parameter ACOUSTIC_FORCING in Par_file'


  !--------------------------------------------------------------------
  !
  ! receivers
  !
  !--------------------------------------------------------------------

  ! read receiver line parameters
  call read_value_integer_p(seismotype, 'solver.seismotype')
  if (err_occurred() /= 0) stop 'error reading parameter seismotype in Par_file'

  call read_value_integer_p(subsamp_seismos, 'solver.subsamp_seismos')
  if (err_occurred() /= 0) stop 'error reading parameter subsamp_seismos in Par_file'

  call read_value_logical_p(USE_TRICK_FOR_BETTER_PRESSURE, 'solver.USE_TRICK_FOR_BETTER_PRESSURE')
  if (err_occurred() /= 0) stop 'error reading parameter USE_TRICK_FOR_BETTER_PRESSURE in Par_file'

  call read_value_integer_p(NSTEP_BETWEEN_OUTPUT_SEISMOS, 'solver.NSTEP_BETWEEN_OUTPUT_SEISMOS')
  if (err_occurred() /= 0) stop 'error reading parameter NSTEP_BETWEEN_OUTPUT_SEISMOS in Par_file'

  call read_value_logical_p(COMPUTE_INTEGRATED_ENERGY_FIELD, 'solver.COMPUTE_INTEGRATED_ENERGY_FIELD')
  if (err_occurred() /= 0) stop 'error reading parameter COMPUTE_INTEGRATED_ENERGY_FIELD in Par_file'

  call read_value_double_precision_p(USER_T0, 'solver.USER_T0')
  if (err_occurred() /= 0) stop 'error reading parameter USER_T0 in Par_file'

  call read_value_logical_p(save_ASCII_seismograms, 'solver.save_ASCII_seismograms')
  if (err_occurred() /= 0) stop 'error reading parameter save_ASCII_seismograms in Par_file'

  call read_value_logical_p(save_binary_seismograms_single, 'solver.save_binary_seismograms_single')
  if (err_occurred() /= 0) stop 'error reading parameter save_binary_seismograms_single in Par_file'

  call read_value_logical_p(save_binary_seismograms_double, 'solver.save_binary_seismograms_double')
  if (err_occurred() /= 0) stop 'error reading parameter save_binary_seismograms_double in Par_file'

  call read_value_logical_p(SU_FORMAT, 'solver.SU_FORMAT')
  if (err_occurred() /= 0) stop 'error reading parameter SU_FORMAT in Par_file'

  call read_value_logical_p(use_existing_STATIONS, 'solver.use_existing_STATIONS')
  if (err_occurred() /= 0) stop 'error reading parameter use_existing_STATIONS in Par_file'

  call read_value_integer_p(nreceiversets, 'solver.nreceiversets')
  if (err_occurred() /= 0) stop 'error reading parameter nreceiversets in Par_file'

  call read_value_double_precision_p(anglerec, 'solver.anglerec')
  if (err_occurred() /= 0) stop 'error reading parameter anglerec in Par_file'

  call read_value_logical_p(rec_normal_to_surface, 'solver.rec_normal_to_surface')
  if (err_occurred() /= 0) stop 'error reading parameter rec_normal_to_surface in Par_file'

  ! reads in receiver sets
  if (nreceiversets < 1) stop 'number of receiver lines must be greater than 1'

  ! allocate receiver line arrays
  allocate(nrec(nreceiversets))
  allocate(xdeb(nreceiversets))
  allocate(zdeb(nreceiversets))
  allocate(xfin(nreceiversets))
  allocate(zfin(nreceiversets))
  allocate(record_at_surface_same_vertical(nreceiversets),stat=ios)
  if (ios /= 0 ) stop 'error allocating receiver lines'

  ! loop on all the receiver lines
  do ireceiverlines = 1,nreceiversets
    call read_value_integer_next_p(nrec(ireceiverlines),'solver.nrec')
    if (err_occurred() /= 0) stop 'error reading parameter nrec in Par_file'

    call read_value_double_prec_next_p(xdeb(ireceiverlines),'solver.xdeb')
    if (err_occurred() /= 0) stop 'error reading parameter xdeb in Par_file'

    call read_value_double_prec_next_p(zdeb(ireceiverlines),'solver.zdeb')
    if (err_occurred() /= 0) stop 'error reading parameter zdeb in Par_file'

    call read_value_double_prec_next_p(xfin(ireceiverlines),'solver.xfin')
    if (err_occurred() /= 0) stop 'error reading parameter xfin in Par_file'

    call read_value_double_prec_next_p(zfin(ireceiverlines),'solver.zfin')
    if (err_occurred() /= 0) stop 'error reading parameter zfin in Par_file'

    call read_value_logical_next_p(record_at_surface_same_vertical(ireceiverlines),'solver.record_at_surface_same_vertical')
    if (err_occurred() /= 0) stop 'error reading parameter record_at_surface_same_vertical in Par_file'

    if (read_external_mesh .and. record_at_surface_same_vertical(ireceiverlines)) then
      stop 'Cannot use record_at_surface_same_vertical with external meshes!'
    endif
  enddo

  !--------------------------------------------------------------------
  !
  ! adjoint kernel outputs
  !
  !--------------------------------------------------------------------

  call read_value_logical_p(save_ASCII_kernels, 'solver.save_ASCII_kernels')
  if (err_occurred() /= 0) stop 'error reading parameter save_ASCII_kernels in Par_file'

  !--------------------------------------------------------------------
  !
  ! boundary conditions
  !
  !--------------------------------------------------------------------

  call read_value_logical_p(PML_BOUNDARY_CONDITIONS, 'solver.PML_BOUNDARY_CONDITIONS')
  if (err_occurred() /= 0) stop 'error reading parameter PML_BOUNDARY_CONDITIONS in Par_file'

  call read_value_integer_p(NELEM_PML_THICKNESS, 'solver.NELEM_PML_THICKNESS')
  if (err_occurred() /= 0) stop 'error reading parameter NELEM_PML_THICKNESS in Par_file'

  call read_value_logical_p(ROTATE_PML_ACTIVATE, 'solver.ROTATE_PML_ACTIVATE')
  if (err_occurred() /= 0) stop 'error reading parameter ROTATE_PML_ACTIVATE in Par_file'

  call read_value_double_precision_p(ROTATE_PML_ANGLE, 'solver.ROTATE_PML_ANGLE')
  if (err_occurred() /= 0) stop 'error reading parameter ROTATE_PML_ANGLE in Par_file'

  ! boolean defining whether to use any absorbing boundaries
  call read_value_logical_p(STACEY_ABSORBING_CONDITIONS, 'solver.STACEY_ABSORBING_CONDITIONS')
  if (err_occurred() /= 0) stop 'error reading parameter STACEY_ABSORBING_CONDITIONS in Par_file'

  call read_value_logical_p(ADD_PERIODIC_CONDITIONS, 'solver.ADD_PERIODIC_CONDITIONS')
  if (err_occurred() /= 0) stop 'error reading parameter ADD_PERIODIC_CONDITIONS in Par_file'

  call read_value_double_precision_p(PERIODIC_HORIZ_DIST, 'solver.PERIODIC_HORIZ_DIST')
  if (err_occurred() /= 0) stop 'error reading parameter PERIODIC_HORIZ_DIST in Par_file'

  !--------------------------------------------------------------------
  !
  ! velocity and density models
  !
  !--------------------------------------------------------------------

  ! read the different material materials (i.e. the number of models)
  call read_value_integer_p(nbmodels, 'mesher.nbmodels')
  if (err_occurred() /= 0) stop 'error reading parameter nbmodels in Par_file'

  ! reads in material definitions
  if (nbmodels <= 0) stop 'Non-positive number of materials not allowed!'

  allocate(icodemat(nbmodels))
  allocate(cp(nbmodels))
  allocate(cs(nbmodels))
  allocate(aniso3(nbmodels))
  allocate(aniso4(nbmodels))
  allocate(aniso5(nbmodels))
  allocate(aniso6(nbmodels))
  allocate(aniso7(nbmodels))
  allocate(aniso8(nbmodels))
  allocate(aniso9(nbmodels))
  allocate(aniso10(nbmodels))
  allocate(aniso11(nbmodels))
  allocate(aniso12(nbmodels))
  allocate(QKappa(nbmodels))
  allocate(Qmu(nbmodels))
  allocate(rho_s(nbmodels))
  allocate(rho_f(nbmodels))
  allocate(phi(nbmodels))
  allocate(tortuosity(nbmodels))
  allocate(permxx(nbmodels))
  allocate(permxz(nbmodels))
  allocate(permzz(nbmodels))
  allocate(kappa_s(nbmodels))
  allocate(kappa_f(nbmodels))
  allocate(kappa_fr(nbmodels))
  allocate(eta_f(nbmodels))
  allocate(mu_fr(nbmodels))

  call read_materials(AXISYM,nbmodels,icodemat,cp,cs, &
                      aniso3,aniso4,aniso5,aniso6,aniso7,aniso8,aniso9,aniso10,aniso11,aniso12, &
                      QKappa,Qmu,rho_s,rho_f,phi,tortuosity, &
                      permxx,permxz,permzz,kappa_s,kappa_f,kappa_fr, &
                      eta_f,mu_fr)

  call read_value_string_p(TOMOGRAPHY_FILE, 'solver.TOMOGRAPHY_FILE')
  if (err_occurred() /= 0) stop 'error reading parameter TOMOGRAPHY_FILE in Par_file'

  ! boolean defining whether internal or external mesh
  call read_value_logical_p(read_external_mesh, 'mesher.read_external_mesh')
  if (err_occurred() /= 0) stop 'error reading parameter read_external_mesh in Par_file'

  !--------------------------------------------------------------------
  !
  ! parameters external / internal meshing
  !
  !--------------------------------------------------------------------


  !-----------------
  ! external mesh parameters

  if (read_external_mesh) then

    ! read info about external mesh
    call read_value_string_p(mesh_file, 'mesher.mesh_file')
    if (err_occurred() /= 0) stop 'error reading parameter mesh_file in Par_file'

    call read_value_string_p(nodes_coords_file, 'mesher.nodes_coords_file')
    if (err_occurred() /= 0) stop 'error reading parameter nodes_coords_file in Par_file'

    call read_value_string_p(materials_file, 'mesher.materials_file')
    if (err_occurred() /= 0) stop 'error reading parameter materials_file in Par_file'

    call read_value_string_p(free_surface_file, 'mesher.free_surface_file')
    if (err_occurred() /= 0) stop 'error reading parameter free_surface_file in Par_file'

    call read_value_string_p(axial_elements_file, 'mesher.axial_elements_file')
    if (err_occurred() /= 0) stop 'error reading parameter axial_elements_file in Par_file'

    call read_value_string_p(absorbing_surface_file, 'mesher.absorbing_surface_file')
    if (err_occurred() /= 0) stop 'error reading parameter absorbing_surface_file in Par_file'

    call read_value_string_p(acoustic_forcing_surface_file, 'mesher.acoustic_forcing_surface_file')
    if (err_occurred() /= 0) stop 'error reading parameter acoustic_forcing_surface_file in Par_file'

    call read_value_string_p(CPML_element_file, 'mesher.CPML_element_file')
    if (err_occurred() /= 0) stop 'error reading parameter CPML_element_file in Par_file'

    call read_value_string_p(tangential_detection_curve_file, 'mesher.tangential_detection_curve_file')
    if (err_occurred() /= 0) stop 'error reading parameter tangential_detection_curve_file in Par_file'

  else

    !-----------------
    ! internal mesh parameters

    ! interfaces file
    call read_value_string_p(interfacesfile, 'mesher.interfacesfile')
    if (err_occurred() /= 0) stop 'error reading parameter interfacesfile in Par_file'

    ! read grid parameters
    call read_value_double_precision_p(xmin_param, 'mesher.xmin')
    if (err_occurred() /= 0) stop 'error reading parameter xmin_param in Par_file'

    call read_value_double_precision_p(xmax_param, 'mesher.xmax')
    if (err_occurred() /= 0) stop 'error reading parameter xmax_param in Par_file'

    call read_value_integer_p(nx_param, 'mesher.nx')
    if (err_occurred() /= 0) stop 'error reading parameter nx in Par_file'

    ! read absorbing boundary parameters
    call read_value_logical_p(absorbbottom, 'solver.absorbbottom')
    if (err_occurred() /= 0) stop 'error reading parameter absorbbottom in Par_file'

    call read_value_logical_p(absorbright, 'solver.absorbright')
    if (err_occurred() /= 0) stop 'error reading parameter absorbright in Par_file'

    call read_value_logical_p(absorbtop, 'solver.absorbtop')
    if (err_occurred() /= 0) stop 'error reading parameter absorbtop in Par_file'

    call read_value_logical_p(absorbleft, 'solver.absorbleft')
    if (err_occurred() /= 0) stop 'error reading parameter absorbleft in Par_file'

    ! note: if internal mesh, then regions will be read in by read_regions (from meshfem2D)

  endif


  !--------------------------------------------------------------------
  !
  ! display parameters
  !
  !--------------------------------------------------------------------

  ! read display parameters
  call read_value_integer_p(NSTEP_BETWEEN_OUTPUT_INFO, 'solver.NSTEP_BETWEEN_OUTPUT_INFO')
  if (err_occurred() /= 0) stop 'error reading parameter NSTEP_BETWEEN_OUTPUT_INFO in Par_file'

  call read_value_logical_p(output_grid_Gnuplot, 'solver.output_grid_Gnuplot')
  if (err_occurred() /= 0) stop 'error reading parameter output_grid_Gnuplot in Par_file'

  call read_value_logical_p(output_grid_ASCII, 'solver.output_grid_ASCII')
  if (err_occurred() /= 0) stop 'error reading parameter output_grid_ASCII in Par_file'

  call read_value_logical_p(output_energy, 'solver.output_energy')
  if (err_occurred() /= 0) stop 'error reading parameter output_energy in Par_file'

  !--------------------------------------------------------------------
  !
  ! movies/images/snapshots
  !
  !--------------------------------------------------------------------

  call read_value_integer_p(NSTEP_BETWEEN_OUTPUT_IMAGES, 'solver.NSTEP_BETWEEN_OUTPUT_IMAGES')
  if (err_occurred() /= 0) stop 'error reading parameter NSTEP_BETWEEN_OUTPUT_IMAGES in Par_file'

  call read_value_integer_p(NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS, 'solver.NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS')
  if (err_occurred() /= 0) stop 'error reading parameter NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS in Par_file'

  call read_value_double_precision_p(cutsnaps, 'solver.cutsnaps')
  if (err_occurred() /= 0) stop 'error reading parameter cutsnaps in Par_file'

  ! jpeg images
  call read_value_logical_p(output_color_image, 'solver.output_color_image')
  if (err_occurred() /= 0) stop 'error reading parameter output_color_image in Par_file'

  call read_value_integer_p(imagetype_JPEG, 'solver.imagetype_JPEG')
  if (err_occurred() /= 0) stop 'error reading parameter imagetype_JPEG in Par_file'

  call read_value_double_precision_p(factor_subsample_image, 'solver.factor_subsample_image')
  if (err_occurred() /= 0) stop 'error reading parameter factor_subsample_image in Par_file'

  call read_value_logical_p(USE_CONSTANT_MAX_AMPLITUDE, 'solver.USE_CONSTANT_MAX_AMPLITUDE')
  if (err_occurred() /= 0) stop 'error reading parameter USE_CONSTANT_MAX_AMPLITUDE in Par_file'

  call read_value_double_precision_p(CONSTANT_MAX_AMPLITUDE_TO_USE, 'solver.CONSTANT_MAX_AMPLITUDE_TO_USE')
  if (err_occurred() /= 0) stop 'error reading parameter CONSTANT_MAX_AMPLITUDE_TO_USE in Par_file'

  call read_value_double_precision_p(POWER_DISPLAY_COLOR, 'solver.POWER_DISPLAY_COLOR')
  if (err_occurred() /= 0) stop 'error reading parameter POWER_DISPLAY_COLOR in Par_file'

  call read_value_logical_p(DRAW_SOURCES_AND_RECEIVERS, 'solver.DRAW_SOURCES_AND_RECEIVERS')
  if (err_occurred() /= 0) stop 'error reading parameter DRAW_SOURCES_AND_RECEIVERS in Par_file'

  call read_value_logical_p(DRAW_WATER_IN_BLUE, 'solver.DRAW_WATER_IN_BLUE')
  if (err_occurred() /= 0) stop 'error reading parameter DRAW_WATER_IN_BLUE in Par_file'

  call read_value_logical_p(USE_SNAPSHOT_NUMBER_IN_FILENAME, 'solver.USE_SNAPSHOT_NUMBER_IN_FILENAME')
  if (err_occurred() /= 0) stop 'error reading parameter USE_SNAPSHOT_NUMBER_IN_FILENAME in Par_file'

  ! postscript files
  call read_value_logical_p(output_postscript_snapshot, 'solver.output_postscript_snapshot')
  if (err_occurred() /= 0) stop 'error reading parameter output_postscript_snapshot in Par_file'

  call read_value_integer_p(imagetype_postscript, 'solver.imagetype_postscript')
  if (err_occurred() /= 0) stop 'error reading parameter imagetype_postscript in Par_file'

  call read_value_logical_p(meshvect, 'solver.meshvect')
  if (err_occurred() /= 0) stop 'error reading parameter meshvect in Par_file'

  call read_value_logical_p(modelvect, 'solver.modelvect')
  if (err_occurred() /= 0) stop 'error reading parameter modelvect in Par_file'

  call read_value_logical_p(boundvect, 'solver.boundvect')
  if (err_occurred() /= 0) stop 'error reading parameter boundvect in Par_file'

  call read_value_logical_p(interpol, 'solver.interpol')
  if (err_occurred() /= 0) stop 'error reading parameter interpol in Par_file'

  call read_value_integer_p(pointsdisp, 'solver.pointsdisp')
  if (err_occurred() /= 0) stop 'error reading parameter pointsdisp in Par_file'

  call read_value_integer_p(subsamp_postscript, 'solver.subsamp_postscript')
  if (err_occurred() /= 0) stop 'error reading parameter subsamp_postscript in Par_file'

  call read_value_double_precision_p(sizemax_arrows, 'solver.sizemax_arrows')
  if (err_occurred() /= 0) stop 'error reading parameter sizemax_arrows in Par_file'

  call read_value_logical_p(US_LETTER, 'solver.US_LETTER')
  if (err_occurred() /= 0) stop 'error reading parameter US_LETTER in Par_file'

  ! wavefield dumps
  call read_value_logical_p(output_wavefield_dumps, 'solver.output_wavefield_dumps')
  if (err_occurred() /= 0) stop 'error reading parameter output_wavefield_dumps in Par_file'

  call read_value_integer_p(imagetype_wavefield_dumps, 'solver.imagetype_wavefield_dumps')
  if (err_occurred() /= 0) stop 'error reading parameter imagetype_wavefield_dumps in Par_file'

  call read_value_logical_p(use_binary_for_wavefield_dumps, 'solver.use_binary_for_wavefield_dumps')
  if (err_occurred() /= 0) stop 'error reading parameter use_binary_for_wavefield_dumps in Par_file'


  ! checks input parameters
  call check_parameters()

  ! boundary conditions
  if (add_Bielak_conditions .and. .not. STACEY_ABSORBING_CONDITIONS) &
    stop 'need STACEY_ABSORBING_CONDITIONS set to .true. in order to use add_Bielak_conditions'

  ! CPML and Stacey are mutually exclusive
  if (STACEY_ABSORBING_CONDITIONS .and. PML_BOUNDARY_CONDITIONS) &
    stop 'STACEY_ABSORBING_CONDITIONS and PML_BOUNDARY_CONDITIONS are mutually exclusive but are both set to .true.'
  ! we also set in subroutine prepare_timerun_read to make sure that STACEY_ABSORBING_CONDITIONS = .false. when
  ! PML_BOUNDARY_CONDITIONS is used.

  ! solve the conflict in value of PML_BOUNDARY_CONDITIONS and STACEY_ABSORBING_CONDITIONS
  if (PML_BOUNDARY_CONDITIONS) STACEY_ABSORBING_CONDITIONS = .true.

  ! initializes flags for absorbing boundaries
  if (.not. STACEY_ABSORBING_CONDITIONS) then
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

  end subroutine read_parameter_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine check_parameters()

  use parameter_file_par

  implicit none

  ! checks partitioning
  if (NPROC <= 0) then
     print *, 'Number of processes (NPROC) must be greater than or equal to one.'
     stop
  endif

#ifndef USE_MPI
  if (NPROC > 1) then
     print *, 'Number of processes (NPROC) must be equal to one when not using MPI.'
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

  ! simulation parameters
  if (SIMULATION_TYPE /= 1 .and. SIMULATION_TYPE /= 3) &
    stop 'SIMULATION_TYPE can only be set to 1 or 3 in the Par_file; exiting'

  if (N_SLS < 2) &
    stop 'must have N_SLS >= 2 even if attenuation if off because it is used to assign some arrays'

  if (ngnod /= 4 .and. ngnod /= 9) &
    stop 'ngnod should be either 4 or 9!'

  if (seismotype < 1 .or. seismotype > 6) &
    stop 'seismotype should be 1(=displ), 2(=veloc), 3(=accel), 4(=pressure), 5(=curl of displ) or 6(=the fluid potential)'

  if (USE_TRICK_FOR_BETTER_PRESSURE .and. seismotype /= 4) &
    stop 'USE_TRICK_FOR_BETTER_PRESSURE : seismograms must record pressure'

  if (subsamp_seismos < 1) &
    stop 'error: subsamp_seismos must be >= 1'

  if (output_color_image .and. USE_CONSTANT_MAX_AMPLITUDE .and. CONSTANT_MAX_AMPLITUDE_TO_USE < 0.d0) &
    stop 'CONSTANT_MAX_AMPLITUDE_TO_USE must be strictly positive'

  if (force_normal_to_surface .or. rec_normal_to_surface) then
    if (.not. read_external_mesh) &
      stop 'Error read_external_mesh must be set to .true. for force_normal_to_surface or rec_normal_to_surface &
            &to use external tangential_dectection_curve_file'
  endif

  ! checks model
  select case (MODEL)
  case ('default','ascii','binary','external','gll','binary_voigt')
    print * ! do nothing
  case default
    stop 'Bad value: MODEL'
  end select

  ! checks model
  select case (SAVE_MODEL)
  case ('default','ascii','binary','external','gll')
    print * ! do nothing
  case default
    stop 'Bad value: SAVE_MODEL'
  end select

  end subroutine check_parameters

