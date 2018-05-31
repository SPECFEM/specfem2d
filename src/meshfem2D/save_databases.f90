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

  subroutine save_databases()

! generates the databases for the solver

  use constants, only: IMAIN,IOUT,MAX_STRING_LEN,OUTPUT_FILES
  use part_unstruct_par, only: nspec,iproc
  use shared_parameters, only: NPROC

  implicit none

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: prname

  do iproc = 0, NPROC-1

    ! opens Database file
    write(prname, "(a,i5.5,a)") './'//trim(OUTPUT_FILES)//'Database',iproc,'.bin'
    open(unit=IOUT,file=trim(prname),status='unknown',action='write',form='unformatted',iostat=ier)
    if (ier /= 0 ) call stop_the_code('Error saving databases; check that directory OUTPUT_FILES exists')

    ! saves header infos and simulation setup
    call save_databases_init()

    ! sources setup
    call save_databases_sources()

    ! mesh node coordinates and mesh properties
    call save_databases_coorg_elem()

    ! attenuation
    call save_databases_attenuation()

    ! material properties
    call save_databases_materials()

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
    close(IOUT)

    ! user output
    if (iproc == 0) write(IMAIN,*)
    write(IMAIN,*) 'slice ',iproc,' has number of spectral elements =',nspec
    call flush_IMAIN()

  enddo

  end subroutine save_databases


!-------------------------------------------------------------------------------

  subroutine save_databases_init()

  use constants, only: IOUT
  use shared_parameters
  use part_unstruct_par, only: iproc,nspec,npgeo,region_pml_external_mesh

  implicit none

  ! prepare total numbers
  ! only counts number of nodes (npgeo) for this partition (iproc)
  call write_glob2loc_nodes_database(IOUT, iproc, npgeo, 1)

  ! only counts number of elements (nspec) for this partition (iproc)
  !   DK DK add support for using PML in MPI mode with external mesh
  !   call write_partition_database(IOUT, iproc, nspec, num_material, ngnod, 1)
  call write_partition_database(IOUT, iproc, nspec, num_material, region_pml_external_mesh, ngnod, 1)

  ! file output
  ! note: older version used ascii-output and annotated all parameters
  !       since this is now in binary format, we remove the text strings as they were not used any further
  ! header infos and simulation setup
  ! '# Database for SPECFEM2D'

  ! 'Title of the simulation'
  write(IOUT) title

  ! 'Type of simulation'
  write(IOUT) NOISE_TOMOGRAPHY, UNDO_ATTENUATION_AND_OR_PML

  ! 'nspec'
  write(IOUT) nspec

  ! 'npgeo NPROC'
  write(IOUT) npgeo,NPROC

  ! 'output_grid_Gnuplot interpol'
  write(IOUT) output_grid_Gnuplot,interpol

  ! 'NSTEP_BETWEEN_OUTPUT_INFO'
  write(IOUT) NSTEP_BETWEEN_OUTPUT_INFO

  ! 'NSTEP_BETWEEN_OUTPUT_SEISMOS'
  write(IOUT) NSTEP_BETWEEN_OUTPUT_SEISMOS

  ! 'NSTEP_BETWEEN_OUTPUT_IMAGES'
  write(IOUT) NSTEP_BETWEEN_OUTPUT_IMAGES

  ! 'PML_BOUNDARY_CONDITIONS'
  write(IOUT) PML_BOUNDARY_CONDITIONS

  ! 'ROTATE_PML_ACTIVATE'
  write(IOUT) ROTATE_PML_ACTIVATE

  ! 'ROTATE_PML_ANGLE'
  write(IOUT) ROTATE_PML_ANGLE

  ! 'K_MIN_PML'
  write(IOUT) K_MIN_PML

  ! 'K_MAX_PML'
  write(IOUT) K_MAX_PML

  ! 'damping_change_factor_acoustic'
  write(IOUT) damping_change_factor_acoustic

  ! 'damping_change_factor_elastic'
  write(IOUT) damping_change_factor_elastic

  ! 'PML_PARAMETER_ADJUSTMENT'
  write(IOUT) PML_PARAMETER_ADJUSTMENT

  ! 'read_external_mesh'
  write(IOUT) read_external_mesh

  ! 'NELEM_PML_THICKNESS'
  write(IOUT) NELEM_PML_THICKNESS

  ! 'subsamp_seismos imagetype_JPEG imagetype_wavefield_dumps'
  write(IOUT) subsamp_seismos,imagetype_JPEG,imagetype_wavefield_dumps

  ! 'output_postscript_snapshot output_color_image'
  write(IOUT) output_postscript_snapshot,output_color_image

  ! 'meshvect modelvect boundvect cutsnaps subsamp_postscript sizemax_arrows'
  write(IOUT) meshvect,modelvect,boundvect,cutsnaps,subsamp_postscript,sizemax_arrows

  ! 'anglerec'
  write(IOUT) anglerec

  ! 'initialfield'
  write(IOUT) initialfield

  ! 'add_Bielak_conditions_bottom','add_Bielak_conditions_right','add_Bielak_conditions_top','add_Bielak_conditions_left'
  write(IOUT) add_Bielak_conditions_bottom,add_Bielak_conditions_right,add_Bielak_conditions_top,add_Bielak_conditions_left

  ! 'seismotype imagetype_postscript'
  write(IOUT) seismotype,imagetype_postscript

  ! 'MODEL'
  write(IOUT) MODEL

  ! 'SAVE_MODEL'
  write(IOUT) SAVE_MODEL

  ! 'TOMOGRAPHY_FILE'
  write(IOUT) TOMOGRAPHY_FILE

  ! 'output_grid_ASCII OUTPUT_ENERGY NTSTEP_BETWEEN_OUTPUT_ENERGY output_wavefield_dumps'
  write(IOUT) output_grid_ASCII,OUTPUT_ENERGY,NTSTEP_BETWEEN_OUTPUT_ENERGY,output_wavefield_dumps

  ! 'use_binary_for_wavefield_dumps'
  write(IOUT) use_binary_for_wavefield_dumps

  ! 'ATTENUATION_VISCOELASTIC ATTENUATION_PORO_FLUID_PART ATTENUATION_VISCOACOUSTIC'
  write(IOUT) ATTENUATION_VISCOELASTIC,ATTENUATION_PORO_FLUID_PART,ATTENUATION_VISCOACOUSTIC

  ! 'USE_SOLVOPT'
  write(IOUT) USE_SOLVOPT

  ! 'save_ASCII_seismograms'
  write(IOUT) save_ASCII_seismograms

  ! 'save_binary_seismograms_single save_binary_seismograms_double'
  write(IOUT) save_binary_seismograms_single,save_binary_seismograms_double

  ! 'USE_TRICK_FOR_BETTER_PRESSURE'
  write(IOUT) USE_TRICK_FOR_BETTER_PRESSURE

  ! 'COMPUTE_INTEGRATED_ENERGY_FIELD'
  write(IOUT) COMPUTE_INTEGRATED_ENERGY_FIELD

  ! 'save_ASCII_kernels'
  write(IOUT) save_ASCII_kernels

  ! 'NSTEP_BETWEEN_COMPUTE_KERNELS'
  write(IOUT) NSTEP_BETWEEN_COMPUTE_KERNELS

  ! 'NO_BACKWARD_RECONSTRUCTION'
  write(IOUT) NO_BACKWARD_RECONSTRUCTION

  ! 'DRAW_SOURCES_AND_RECEIVERS'
  write(IOUT) DRAW_SOURCES_AND_RECEIVERS

  ! 'Q0_poroelastic freq0_poroelastic'
  write(IOUT) Q0_poroelastic,freq0_poroelastic

  ! 'Axisymmetric (2.5D, .true.) or Cartesian planar (2D; .false.) simulation'
  write(IOUT) AXISYM

  ! 'P_SV'
  write(IOUT) P_SV

  ! 'factor_subsample_image'
  write(IOUT) factor_subsample_image

  ! 'USE_CONSTANT_MAX_AMPLITUDE'
  write(IOUT) USE_CONSTANT_MAX_AMPLITUDE

  ! 'CONSTANT_MAX_AMPLITUDE_TO_USE'
  write(IOUT) CONSTANT_MAX_AMPLITUDE_TO_USE

  ! 'USE_SNAPSHOT_NUMBER_IN_FILENAME'
  write(IOUT) USE_SNAPSHOT_NUMBER_IN_FILENAME

  ! 'DRAW_WATER_IN_BLUE'
  write(IOUT) DRAW_WATER_IN_BLUE

  ! 'US_LETTER'
  write(IOUT) US_LETTER

  ! 'POWER_DISPLAY_COLOR'
  write(IOUT) POWER_DISPLAY_COLOR

  ! 'SU_FORMAT'
  write(IOUT) SU_FORMAT

  ! 'USER_T0'
  write(IOUT) USER_T0

  ! 'time_stepping_scheme'
  write(IOUT) time_stepping_scheme

  ! 'ADD_PERIODIC_CONDITIONS'
  write(IOUT) ADD_PERIODIC_CONDITIONS

  ! 'PERIODIC_HORIZ_DIST'
  write(IOUT) PERIODIC_HORIZ_DIST

  ! 'GPU_MODE'
  write(IOUT) GPU_MODE

  ! 'setup_with_binary_database"
  write(IOUT) setup_with_binary_database

  ! 'NSTEP DT'
  write(IOUT) NSTEP,DT

  ! 'NT_DUMP_ATTENUATION'
  write(IOUT) NT_DUMP_ATTENUATION

  ! 'ACOUSTIC_FORCING'
  write(IOUT) ACOUSTIC_FORCING

  ! 'NUMBER_OF_SIMULTANEOUS_RUNS'
  write(IOUT) NUMBER_OF_SIMULTANEOUS_RUNS

  ! 'BROADCAST_SAME_MESH_AND_MODEL'
  write(IOUT) BROADCAST_SAME_MESH_AND_MODEL

  end subroutine save_databases_init

!-------------------------------------------------------------------------------

  subroutine save_databases_sources()

  use constants, only: IOUT
  use source_file_par
  use shared_parameters

  implicit none
  ! local parameters
  integer :: i_source

  ! sources setup
  ! 'NSOURCES'
  write(IOUT) NSOURCES

  do i_source = 1,NSOURCES
    write(IOUT) source_type(i_source),time_function_type(i_source)
    write(IOUT) name_of_source_file(i_source)
    write(IOUT) burst_band_width(i_source)
    write(IOUT) xs(i_source),zs(i_source)
    write(IOUT) f0_source(i_source),tshift_src(i_source)
    write(IOUT) factor(i_source),anglesource(i_source)
    write(IOUT) Mxx(i_source),Mzz(i_source),Mxz(i_source)
  enddo

  end subroutine save_databases_sources


!-------------------------------------------------------------------------------

  subroutine save_databases_coorg_elem()

  use constants, only: IOUT
  use part_unstruct_par
  use shared_parameters

  implicit none

  ! 'Coordinates of macrobloc mesh (coorg):'

  call write_glob2loc_nodes_database(IOUT, iproc, npgeo, 2)

  ! 'numat ngnod nspec pointsdisp plot_lowerleft_corner_only'
  write(IOUT) nbmodels,ngnod,nspec,pointsdisp,plot_lowerleft_corner_only

  ! counts number of absorbing elements
  if (any_abs) then
    call write_abs_merge_database(IOUT, iproc, 1)
  else
    nelemabs_loc = 0
  endif

  ! counts number of acoustic forcing elements
  if (ACOUSTIC_FORCING) then
    call write_acoustic_forcing_merge_database(IOUT, iproc, 1)
  else
    nelemacforcing_loc = 0
  endif

  ! counts number of acoustic surface elements
  call write_surface_database(IOUT, nelem_acoustic_surface, acoustic_surface, nelem_acoustic_surface_loc,iproc, 1)

  ! counts number of coupling edges
  ! fluid-solid
  call write_fluidsolid_edges_database(IOUT, nedges_coupled, edges_coupled, nedges_coupled_loc, iproc, 1)
  ! fluid-poroelastic
  call write_fluidsolid_edges_database(IOUT, nedges_acporo_coupled, edges_acporo_coupled, nedges_acporo_coupled_loc, iproc, 1)
  ! solid-poroelastic
  call write_fluidsolid_edges_database(IOUT, nedges_elporo_coupled, edges_elporo_coupled, nedges_elporo_coupled_loc, iproc, 1)

  ! counts number of axial elements
  call write_axial_elements_database(IOUT, nelem_on_the_axis, ispec_of_axial_elements, &
                                     nelem_on_the_axis_loc, iproc, 1, remove_min_to_start_at_zero)

  ! 'nelemabs nelemacforcing nelem_acoustic_surface num_fluid_solid_edges'
  ! 'num_fluid_poro_edges num_solid_poro_edges'
  ! 'nnodes_tangential_curve nelem_on_the_axis'
  write(IOUT) nelemabs_loc,nelemacforcing_loc,nelem_acoustic_surface_loc, &
              nedges_coupled_loc,nedges_acporo_coupled_loc, &
              nedges_elporo_coupled_loc,nnodes_tangential_curve, &
              nelem_on_the_axis_loc

  end subroutine save_databases_coorg_elem

!-------------------------------------------------------------------------------

  subroutine save_databases_attenuation()

  use constants, only: IOUT
  use part_unstruct_par
  use shared_parameters

  implicit none

  ! attenuation setting
  ! 'attenuation'
  write(IOUT) N_SLS, ATTENUATION_f0_REFERENCE, READ_VELOCITIES_AT_f0

  end subroutine save_databases_attenuation

!-------------------------------------------------------------------------------

  subroutine save_databases_materials()

  use constants, only: IOUT,ISOTROPIC_MATERIAL,ANISOTROPIC_MATERIAL,POROELASTIC_MATERIAL
  use part_unstruct_par
  use shared_parameters

  implicit none

  ! local parameters
  double precision :: val0,val1,val2,val3,val4,val5,val6,val7,val8,val9,val10,val11,val12
  integer :: i,indic

  ! material set header
  ! 'Material sets (num 1 rho vp vs 0 0 QKappa Qmu 0 0 0 0 0 0) or '
  ! '(num 2 rho c11 c13 c15 c33 c35 c55 c12 c23 c25 0 0 0) or '
  ! '(num 3 rhos rhof phi c k_xx k_xz k_zz Ks Kf Kfr etaf mufr Qmu)'

  do i = 1,nbmodels
    ! material type
    indic = icodemat(i)

    if (indic == ISOTROPIC_MATERIAL) then
      ! isotropic
      val0 = rho_s_read(i)
      val1 = cp(i)
      val2 = cs(i)
      val3 = comp_g(i)
      val4 = 0.d0
      val5 = QKappa(i)
      val6 = Qmu(i)
      val7 = 0.d0
      val8 = 0.d0
      val9 = 0.d0
      val10 = 0.d0
      val11 = 0.d0
      val12 = 0.d0
      ! old
      !write(IOUT) i,icodemat(i),rho_s_read(i),cp(i),cs(i),0,0,QKappa(i),Qmu(i),0,0,0,0,0,0

    else if (indic == ANISOTROPIC_MATERIAL) then
      ! anisotropic
      val0 = rho_s_read(i)
      val1 = aniso3(i)
      val2 = aniso4(i)
      val3 = aniso5(i)
      val4 = aniso6(i)
      val5 = aniso7(i)
      val6 = aniso8(i)
      val7 = aniso9(i)
      val8 = aniso10(i)
      val9 = aniso11(i)
      val10 = aniso12(i)
      val11 = 0.d0
      val12 = 0.d0
      ! old
      !write(IOUT) i,icodemat(i),rho_s_read(i), &
      !            aniso3(i),aniso4(i),aniso5(i),aniso6(i), &
      !            aniso7(i),aniso8(i),aniso9(i),aniso10(i),aniso11(i),aniso12(i),0,0

    else if (indic == POROELASTIC_MATERIAL) then
      ! poro-elastic
      val0 = rho_s_read(i)
      val1 = rho_f_read(i)
      val2 = phi_read(i)
      val3 = tortuosity_read(i)
      val4 = permxx_read(i)
      val5 = permxz_read(i)
      val6 = permzz_read(i)
      val7 = kappa_s_read(i)
      val8 = kappa_f_read(i)
      val9 = kappa_fr_read(i)
      val10 = eta_f_read(i)
      val11 = mu_fr_read(i)
      val12 = Qmu(i)
      ! old
      !write(IOUT) i,icodemat(i),rho_s_read(i),rho_f_read(i),phi_read(i),tortuosity_read(i), &
      !            permxx_read(i),permxz_read(i),permzz_read(i),kappa_s_read(i), &
      !            kappa_f_read(i),kappa_fr_read(i),eta_f_read(i),mu_fr_read(i),Qmu(i)

    else if (indic <= 0) then
      ! external material
      ! The values will be read from an external tomo file
      val0 = rho_s_read(i)
      val1 = cp(i)
      val2 = cs(i)
      val3 = 0.d0
      val4 = 0.d0
      val5 = QKappa(i)
      val6 = Qmu(i)
      val7 = 0.d0
      val8 = 0.d0
      val9 = 0.d0
      val10 = 0.d0
      val11 = 0.d0
      val12 = 0.d0
      ! old
      !write(IOUT) i,icodemat(i),rho_s_read(i),cp(i),cs(i),0,0,QKappa(i),Qmu(i),0,0,0,0,0,0
    else
      ! case should not occur
      call stop_the_code('Unknown material code')
    endif

    ! check format with file src/specfem2D/read_materials.f90
    write(IOUT) i,indic,val0,val1,val2,val3,val4,val5,val6,val7,val8,val9,val10,val11,val12
  enddo

  ! writes out material properties
  ! 'Arrays kmato and knods for each bloc:'

!   DK DK add support for using PML in MPI mode with external mesh
!   call write_partition_database(IOUT, iproc, nspec, num_material, ngnod, 2)
  call write_partition_database(IOUT, iproc, nspec, num_material, region_pml_external_mesh, ngnod, 2)

  end subroutine save_databases_materials

!-------------------------------------------------------------------------------

  subroutine save_databases_interfaces()

  use constants, only: IOUT
  use part_unstruct_par
  use decompose_par
  use shared_parameters, only: NPROC

  implicit none

  ! local parameters
  integer :: idummy
  integer :: max_interfaces
  integer :: my_ninterface

  if (NPROC /= 1) then
    ! counts interfaces
    call write_interfaces_database(IOUT, NPROC, iproc,my_ninterface, my_interfaces, my_nb_interfaces, 1)

    max_interfaces = maxval(my_nb_interfaces)

    ! 'Interfaces:'
    write(IOUT) my_ninterface, max_interfaces

    ! writes out interface infos
    call write_interfaces_database(IOUT, NPROC, iproc,my_ninterface, my_interfaces, my_nb_interfaces, 2)

  else
    ! single partition, no interfaces
    ! dummy
    idummy = 0
    ! 'Interfaces:'
    write(IOUT) idummy, idummy
  endif

  end subroutine save_databases_interfaces

!-------------------------------------------------------------------------------

  subroutine save_databases_absorbing()

  use constants, only: IOUT
  use shared_parameters, only: any_abs
  use part_unstruct_par, only: iproc

  implicit none

  ! 'List of absorbing elements (edge1 edge2 edge3 edge4 type):'
  if (any_abs) then
    ! writes out absorbing boundaries
    call write_abs_merge_database(IOUT, iproc, 2)
  endif

  end subroutine save_databases_absorbing

!-------------------------------------------------------------------------------

  subroutine save_databases_acoustic_forcing()

  use constants, only: IOUT
  use shared_parameters, only: ACOUSTIC_FORCING
  use part_unstruct_par, only: iproc

  implicit none

  ! 'List of acoustic forcing elements (edge1 edge2 edge3 edge4 type):'
  if (ACOUSTIC_FORCING) then
    ! writes out acoustic forcing edges
    call write_acoustic_forcing_merge_database(IOUT, iproc, 2)
  endif

  end subroutine save_databases_acoustic_forcing

!-------------------------------------------------------------------------------

  subroutine save_databases_free_surf()

  use constants, only: IOUT
  use part_unstruct_par

  implicit none

  ! 'List of acoustic free-surface elements:'
  call write_surface_database(IOUT, nelem_acoustic_surface, acoustic_surface, nelem_acoustic_surface_loc,iproc, 2)

  end subroutine save_databases_free_surf

!-------------------------------------------------------------------------------

  subroutine save_databases_coupled()

  use constants, only: IOUT
  use part_unstruct_par

  implicit none

  ! 'List of acoustic elastic coupled edges:'
  call write_fluidsolid_edges_database(IOUT, nedges_coupled, edges_coupled, nedges_coupled_loc, iproc, 2)

  ! 'List of acoustic poroelastic coupled edges:'
  call write_fluidsolid_edges_database(IOUT, nedges_acporo_coupled, edges_acporo_coupled, nedges_acporo_coupled_loc, iproc, 2)

  ! 'List of poroelastic elastic coupled edges:'
  call write_fluidsolid_edges_database(IOUT, nedges_elporo_coupled, edges_elporo_coupled, nedges_elporo_coupled_loc, iproc, 2)

  end subroutine save_databases_coupled

!-------------------------------------------------------------------------------

  subroutine save_databases_tangential()

  use constants, only: IOUT
  use part_unstruct_par
  use shared_parameters

  implicit none

  ! local parameters
  integer :: i

  ! 'List of tangential detection curve nodes:'
  write(IOUT) force_normal_to_surface,rec_normal_to_surface

  if (force_normal_to_surface .or. rec_normal_to_surface) then
    do i = 1, nnodes_tangential_curve
      write(IOUT) nodes_tangential_curve(1,i),nodes_tangential_curve(2,i)
    enddo
  endif

  end subroutine save_databases_tangential

!-------------------------------------------------------------------------------

  subroutine save_databases_axial_elements

  use constants, only: IOUT
  use part_unstruct_par

  implicit none

  ! 'List of axial elements:'
  call write_axial_elements_database(IOUT, nelem_on_the_axis, ispec_of_axial_elements, &
                                     nelem_on_the_axis_loc, iproc, 2, remove_min_to_start_at_zero)

  end subroutine save_databases_axial_elements



