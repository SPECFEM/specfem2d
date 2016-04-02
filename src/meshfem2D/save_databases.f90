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


  subroutine save_databases(nspec,region_pml_external_mesh,remove_min_to_start_at_zero)

! generates the databases for the solver


  use part_unstruct_par,only: nelmnts
  use parameter_file_par,only: NPROC

  implicit none

  include "constants.h"

  integer :: nspec,remove_min_to_start_at_zero
  integer, dimension(nelmnts) :: region_pml_external_mesh

  ! local parameters
  integer :: iproc,i,ier
  integer :: npgeo
  integer :: my_ninterface
  integer :: nedges_coupled_loc
  integer :: nedges_acporo_coupled_loc
  integer :: nedges_elporo_coupled_loc
  character(len=MAX_STRING_LEN) :: prname


  do iproc = 0, NPROC-1

    ! opens Database file
    write(prname, "('./OUTPUT_FILES/Database',i5.5)") iproc
    open(unit=15,file=trim(prname),status='unknown',iostat=ier)
    if (ier /= 0 ) stop 'Error saving databases; check that directory OUTPUT_FILES exists'

    ! saves header infos and simulation setup
    call save_databases_init()

    ! sources setup
    call save_databases_sources()

    ! mesh node coordinates and mesh properties
    call save_databases_coorg_elem()

    ! attenuation
    call save_databases_attenuation()

    ! material properties
    call save_databases_kmato()

    ! saves mpi/partition interface
    call save_databases_interfaces()

    ! absorbing boundary
    call save_databases_absorbing()

    ! acoustic forcing
    call save_databases_acoustic_forcing()

    ! acoustic free surface
    call save_databases_free_surf()

    ! coupled domain edges
    call save_databases_coupled()

    ! tangential detection curve
    call save_databases_tangential()

    ! axial elements
    call save_databases_axial_elements()

    ! closes Database file
    close(15)

    ! user output
    if (iproc == 0) print *,''
    print *,'slice ',iproc,' has number of spectral elements =',nspec

  enddo

  contains

!-------------------------------------------------------------------------------

  subroutine save_databases_init()

  use parameter_file_par

  implicit none

  ! header infos and simulation setup
  write(15,*) '#'
  write(15,*) '# Database for SPECFEM2D'
  write(15,*) '# Dimitri Komatitsch, (c) CNRS, France'
  write(15,*) '#'

  write(15,*) 'Title of the simulation'
  write(15,'(a100)') title

  write(15,*) 'Axisymmetric (2.5D, .true.) or Cartesian planar (2D; .false.) simulation'
  write(15,*) AXISYM

  write(15,*) 'Type of simulation'
  write(15,*) SIMULATION_TYPE, NOISE_TOMOGRAPHY, SAVE_FORWARD, UNDO_ATTENUATION

  ! counts number of nodes (npgeo) for this partition (iproc)
  call write_glob2loc_nodes_database(15, iproc, npgeo, 1)

  ! counts number of elements (nspec) for this partition (iproc)
  !   DK DK add support for using pml in mpi mode with external mesh
  !   call write_partition_database(15, iproc, nspec, num_material, ngnod, 1)
  call write_partition_database(15, iproc, nspec, num_material, region_pml_external_mesh, ngnod, 1)

  write(15,*) 'nspec'
  write(15,*) nspec

  write(15,*) 'npgeo NPROC'
  write(15,*) npgeo,NPROC

  write(15,*) 'output_grid_Gnuplot interpol'
  write(15,*) output_grid_Gnuplot,interpol

  write(15,*) 'NSTEP_BETWEEN_OUTPUT_INFO'
  write(15,*) NSTEP_BETWEEN_OUTPUT_INFO

  write(15,*) 'NSTEP_BETWEEN_OUTPUT_SEISMOS'
  write(15,*) NSTEP_BETWEEN_OUTPUT_SEISMOS

  write(15,*) 'NSTEP_BETWEEN_OUTPUT_IMAGES'
  write(15,*) NSTEP_BETWEEN_OUTPUT_IMAGES

  write(15,*) 'PML_BOUNDARY_CONDITIONS'
  write(15,*) PML_BOUNDARY_CONDITIONS

  write(15,*) 'ROTATE_PML_ACTIVATE'
  write(15,*) ROTATE_PML_ACTIVATE

  write(15,*) 'ROTATE_PML_ANGLE'
  write(15,*) ROTATE_PML_ANGLE

  write(15,*) 'read_external_mesh'
  write(15,*) read_external_mesh

  write(15,*) 'NELEM_PML_THICKNESS'
  write(15,*) NELEM_PML_THICKNESS

  write(15,*) 'NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS'
  write(15,*) NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS

  write(15,*) 'subsamp_seismos imagetype_JPEG imagetype_wavefield_dumps'
  write(15,*) subsamp_seismos,imagetype_JPEG,imagetype_wavefield_dumps

  write(15,*) 'output_postscript_snapshot output_color_image colors numbers'
  write(15,*) output_postscript_snapshot,output_color_image,' 1 0'

  write(15,*) 'meshvect modelvect boundvect cutsnaps subsamp_postscript sizemax_arrows'
  write(15,*) meshvect,modelvect,boundvect,cutsnaps,subsamp_postscript,sizemax_arrows

  write(15,*) 'anglerec'
  write(15,*) anglerec

  write(15,*) 'initialfield'
  write(15,*) initialfield

  write(15,*) 'add_Bielak_conditions_bottom'
  write(15,*) add_Bielak_conditions_bottom

  write(15,*) 'add_Bielak_conditions_right'
  write(15,*) add_Bielak_conditions_right

  write(15,*) 'add_Bielak_conditions_top'
  write(15,*) add_Bielak_conditions_top

  write(15,*) 'add_Bielak_conditions_left'
  write(15,*) add_Bielak_conditions_left

  write(15,*) 'seismotype imagetype_postscript'
  write(15,*) seismotype,imagetype_postscript

  write(15,*) 'MODEL'
  write(15,'(a100)') MODEL

  write(15,*) 'SAVE_MODEL'
  write(15,'(a100)') SAVE_MODEL

  write(15,*) 'TOMOGRAPHY_FILE'
  write(15,'(a100)') TOMOGRAPHY_FILE

  write(15,*) 'output_grid_ASCII output_energy output_wavefield_dumps'
  write(15,*) output_grid_ASCII,output_energy,output_wavefield_dumps

  write(15,*) 'use_binary_for_wavefield_dumps'
  write(15,*) use_binary_for_wavefield_dumps

  write(15,*) 'ATTENUATION_VISCOELASTIC_SOLID ATTENUATION_PORO_FLUID_PART'
  write(15,*) ATTENUATION_VISCOELASTIC_SOLID,ATTENUATION_PORO_FLUID_PART

  write(15,*) 'save_ASCII_seismograms'
  write(15,*) save_ASCII_seismograms

  write(15,*) 'save_binary_seismograms_single save_binary_seismograms_double'
  write(15,*) save_binary_seismograms_single,save_binary_seismograms_double

  write(15,*) 'USE_TRICK_FOR_BETTER_PRESSURE'
  write(15,*) USE_TRICK_FOR_BETTER_PRESSURE

  write(15,*) 'COMPUTE_INTEGRATED_ENERGY_FIELD'
  write(15,*) COMPUTE_INTEGRATED_ENERGY_FIELD

  write(15,*) 'save_ASCII_kernels'
  write(15,*) save_ASCII_kernels

  write(15,*) 'DRAW_SOURCES_AND_RECEIVERS'
  write(15,*) DRAW_SOURCES_AND_RECEIVERS

  write(15,*) 'Q0 freq0'
  write(15,*) Q0,freq0

  write(15,*) 'P_SV'
  write(15,*) P_SV

  write(15,*) 'factor_subsample_image'
  write(15,*) factor_subsample_image

  write(15,*) 'USE_CONSTANT_MAX_AMPLITUDE'
  write(15,*) USE_CONSTANT_MAX_AMPLITUDE

  write(15,*) 'CONSTANT_MAX_AMPLITUDE_TO_USE'
  write(15,*) CONSTANT_MAX_AMPLITUDE_TO_USE

  write(15,*) 'USE_SNAPSHOT_NUMBER_IN_FILENAME'
  write(15,*) USE_SNAPSHOT_NUMBER_IN_FILENAME

  write(15,*) 'DRAW_WATER_IN_BLUE'
  write(15,*) DRAW_WATER_IN_BLUE

  write(15,*) 'US_LETTER'
  write(15,*) US_LETTER

  write(15,*) 'POWER_DISPLAY_COLOR'
  write(15,*) POWER_DISPLAY_COLOR

  write(15,*) 'SU_FORMAT'
  write(15,*) SU_FORMAT

  write(15,*) 'USER_T0'
  write(15,*) USER_T0

  write(15,*) 'time_stepping_scheme'
  write(15,*) time_stepping_scheme

  write(15,*) 'ADD_PERIODIC_CONDITIONS'
  write(15,*) ADD_PERIODIC_CONDITIONS

  write(15,*) 'PERIODIC_HORIZ_DIST'
  write(15,*) PERIODIC_HORIZ_DIST

  write(15,*) 'GPU_MODE'
  write(15,*) GPU_MODE

  write(15,*) 'NSTEP DT'
  write(15,*) NSTEP,DT

  write(15,*) 'NT_DUMP_ATTENUATION'
  write(15,*) NT_DUMP_ATTENUATION

  write(15,*) 'ACOUSTIC_FORCING'
  write(15,*) ACOUSTIC_FORCING

  end subroutine save_databases_init

!-------------------------------------------------------------------------------

  subroutine save_databases_sources()

  use source_file_par
  use parameter_file_par

  implicit none
  ! local parameters
  integer :: i_source

  ! sources setup
  write(15,*) 'NSOURCES'
  write(15,*) NSOURCES

  do i_source = 1,NSOURCES
    write(15,*) 'source', i_source
    write(15,*) source_type(i_source),time_function_type(i_source)
    write(15,'(a100)') name_of_source_file(i_source) ! aXXX: Write wrong character if XXX != 100 !!!
    write(15,*) burst_band_width(i_source), &
                xs(i_source),zs(i_source),f0_source(i_source),tshift_src(i_source), &
                factor(i_source),anglesource(i_source), &
                Mxx(i_source),Mzz(i_source),Mxz(i_source)
  enddo

  end subroutine save_databases_sources


!-------------------------------------------------------------------------------

  subroutine save_databases_coorg_elem()

  use part_unstruct_par
  use parameter_file_par

  implicit none

  write(15,*) 'Coordinates of macrobloc mesh (coorg):'

  call write_glob2loc_nodes_database(15, iproc, npgeo, 2)

  write(15,*) 'numat ngnod nspec pointsdisp plot_lowerleft_corner_only'
  write(15,*) nbmodels,ngnod,nspec,pointsdisp,plot_lowerleft_corner_only

  ! counts number of absorbing elements
  if (STACEY_ABSORBING_CONDITIONS) then
    call write_abs_merge_database(15, iproc, 1)
  else
    nelemabs_loc = 0
  endif

  ! counts number of acoustic forcing elements
  if (ACOUSTIC_FORCING) then
    call write_acoustic_forcing_merge_database(15, iproc, 1)
  else
    nelemacforcing_loc = 0
  endif

  ! counts number of acoustic surface elements
  call write_surface_database(15, nelem_acoustic_surface, acoustic_surface, nelem_acoustic_surface_loc,iproc, 1)

  ! counts number of coupling edges
  ! fluid-solid
  call write_fluidsolid_edges_database(15, nedges_coupled, edges_coupled, nedges_coupled_loc, iproc, 1)
  ! fluid-poroelastic
  call write_fluidsolid_edges_database(15, nedges_acporo_coupled, edges_acporo_coupled, nedges_acporo_coupled_loc, iproc, 1)
  ! solid-poroelastic
  call write_fluidsolid_edges_database(15, nedges_elporo_coupled, edges_elporo_coupled, nedges_elporo_coupled_loc, iproc, 1)

  ! counts number of axial elements
  call write_axial_elements_database(15, nelem_on_the_axis, ispec_of_axial_elements, &
                                     nelem_on_the_axis_loc, iproc, 1, remove_min_to_start_at_zero)

  write(15,*) 'nelemabs nelemacforcing nelem_acoustic_surface num_fluid_solid_edges'
  write(15,*) 'num_fluid_poro_edges num_solid_poro_edges'
  write(15,*) 'nnodes_tangential_curve nelem_on_the_axis'
  write(15,*) nelemabs_loc,nelemacforcing_loc,nelem_acoustic_surface_loc, &
              nedges_coupled_loc,nedges_acporo_coupled_loc,&
              nedges_elporo_coupled_loc,nnodes_tangential_curve, &
              nelem_on_the_axis_loc

  end subroutine save_databases_coorg_elem

!-------------------------------------------------------------------------------

  subroutine save_databases_attenuation()

  use part_unstruct_par
  use parameter_file_par

  implicit none

  ! attenuation setting
  write(15,*) 'attenuation'
  write(15,*) N_SLS, f0_attenuation, READ_VELOCITIES_AT_f0

  end subroutine save_databases_attenuation

!-------------------------------------------------------------------------------

  subroutine save_databases_kmato()

  use part_unstruct_par
  use parameter_file_par

  implicit none

  ! material set header
  write(15,*) 'Material sets (num 1 rho vp vs 0 0 QKappa Qmu 0 0 0 0 0 0) or '
  write(15,*) '(num 2 rho c11 c13 c15 c33 c35 c55 c12 c23 c25 0 0 0) or '
  write(15,*) '(num 3 rhos rhof phi c k_xx k_xz k_zz Ks Kf Kfr etaf mufr Qmu)'

  do i = 1,nbmodels
    if (icodemat(i) == ISOTROPIC_MATERIAL) then
      ! isotropic
      write(15,*) i,icodemat(i),rho_s(i),cp(i),cs(i),0,0,QKappa(i),Qmu(i),0,0,0,0,0,0

    else if (icodemat(i) == ANISOTROPIC_MATERIAL) then
      ! anisotropic
      write(15,*) i,icodemat(i),rho_s(i), &
                  aniso3(i),aniso4(i),aniso5(i),aniso6(i),&
                  aniso7(i),aniso8(i),aniso9(i),aniso10(i),aniso11(i),aniso12(i),0,0

    else if (icodemat(i) == POROELASTIC_MATERIAL) then
      ! poro-elastic
      write(15,*) i,icodemat(i),rho_s(i),rho_f(i),phi(i),tortuosity(i), &
                  permxx(i),permxz(i),permzz(i),kappa_s(i),&
                  kappa_f(i),kappa_fr(i),eta_f(i),mu_fr(i),Qmu(i)

    else if (icodemat(i) <= 0) then
      ! external material
      ! The values will be read from an external tomo file
      write(15,*) i,icodemat(i),rho_s(i),cp(i),cs(i),0,0,QKappa(i),Qmu(i),0,0,0,0,0,0
    else
      ! case should not occur
      stop 'Unknown material code'
    endif
  enddo

  ! writes out material properties
  write(15,*) 'Arrays kmato and knods for each bloc:'

!   DK DK add support for using PML in MPI mode with external mesh
!   call write_partition_database(15, iproc, nspec, num_material, ngnod, 2)
  call write_partition_database(15, iproc, nspec, num_material, region_pml_external_mesh, ngnod, 2)

  end subroutine save_databases_kmato

!-------------------------------------------------------------------------------

  subroutine save_databases_interfaces()

  use part_unstruct_par
  use decompose_par

  implicit none

  if (NPROC /= 1) then
    ! counts interfaces
    call write_interfaces_database(15, NPROC, iproc,my_ninterface, my_interfaces, my_nb_interfaces, 1)

    write(15,*) 'Interfaces:'
    write(15,*) my_ninterface, maxval(my_nb_interfaces)

    ! writes out interface infos
    call write_interfaces_database(15, NPROC, iproc,my_ninterface, my_interfaces, my_nb_interfaces, 2)

  else
    ! single partition, no interfaces
    ! dummy
    write(15,*) 'Interfaces:'
    write(15,*) 0, 0
  endif

  end subroutine save_databases_interfaces

!-------------------------------------------------------------------------------

  subroutine save_databases_absorbing()

  use parameter_file_par,only: STACEY_ABSORBING_CONDITIONS

  implicit none

  write(15,*) 'List of absorbing elements (edge1 edge2 edge3 edge4 type):'
  if (STACEY_ABSORBING_CONDITIONS) then
    ! writes out absorbing boundaries
    call write_abs_merge_database(15, iproc, 2)
  endif

  end subroutine save_databases_absorbing

!-------------------------------------------------------------------------------

  subroutine save_databases_acoustic_forcing()

  use parameter_file_par,only: ACOUSTIC_FORCING

  implicit none

  write(15,*) 'List of acoustic forcing elements (edge1 edge2 edge3 edge4 type):'
  if (ACOUSTIC_FORCING) then
    ! writes out acoustic forcing edges
    call write_acoustic_forcing_merge_database(15, iproc, 2)
  endif

  end subroutine save_databases_acoustic_forcing

!-------------------------------------------------------------------------------

  subroutine save_databases_free_surf()

  use part_unstruct_par

  implicit none

  write(15,*) 'List of acoustic free-surface elements:'
  call write_surface_database(15, nelem_acoustic_surface, acoustic_surface, nelem_acoustic_surface_loc,iproc, 2)

  end subroutine save_databases_free_surf

!-------------------------------------------------------------------------------

  subroutine save_databases_coupled()

  use part_unstruct_par

  implicit none

  write(15,*) 'List of acoustic elastic coupled edges:'
  call write_fluidsolid_edges_database(15, nedges_coupled, edges_coupled, nedges_coupled_loc, iproc, 2)

  write(15,*) 'List of acoustic poroelastic coupled edges:'
  call write_fluidsolid_edges_database(15, nedges_acporo_coupled, edges_acporo_coupled, nedges_acporo_coupled_loc, iproc, 2)

  write(15,*) 'List of poroelastic elastic coupled edges:'
  call write_fluidsolid_edges_database(15, nedges_elporo_coupled, edges_elporo_coupled, nedges_elporo_coupled_loc, iproc, 2)

  end subroutine save_databases_coupled

!-------------------------------------------------------------------------------

  subroutine save_databases_tangential()

  use part_unstruct_par
  use parameter_file_par

  implicit none

  write(15,*) 'List of tangential detection curve nodes:'
  !write(15,*) nnodes_tangential_curve
  write(15,*) force_normal_to_surface,rec_normal_to_surface

  if (force_normal_to_surface .or. rec_normal_to_surface) then
    do i = 1, nnodes_tangential_curve
      write(15,*) nodes_tangential_curve(1,i),nodes_tangential_curve(2,i)
    enddo
  endif

  end subroutine save_databases_tangential

!-------------------------------------------------------------------------------

  subroutine save_databases_axial_elements

  use part_unstruct_par

  implicit none

  write(15,*) 'List of axial elements:'
  call write_axial_elements_database(15, nelem_on_the_axis, ispec_of_axial_elements, &
                                     nelem_on_the_axis_loc, iproc, 2, remove_min_to_start_at_zero)

  end subroutine save_databases_axial_elements

!-------------------------------------------------------------------------------

  end subroutine save_databases

