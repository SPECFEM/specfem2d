
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
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and Inria at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================


  subroutine save_databases(nspec,num_material,region_pml_external_mesh, &
                            my_interfaces,my_nb_interfaces, &
                            nnodes_tangential_curve,nodes_tangential_curve,remove_min_to_start_at_zero)


! generates the databases for the solver

  use part_unstruct
  use parameter_file
  use source_file
  implicit none
  include "constants.h"

  integer :: nspec,remove_min_to_start_at_zero
  integer, dimension(nelmnts) :: num_material
  integer, dimension(nelmnts) :: region_pml_external_mesh

  integer, dimension(0:ninterfaces-1) :: my_interfaces
  integer, dimension(0:ninterfaces-1) :: my_nb_interfaces

  integer ::  nnodes_tangential_curve
  double precision, dimension(2,nnodes_tangential_curve) :: nodes_tangential_curve

  ! local parameters
  integer :: iproc,i_source,i,ios
  integer :: npgeo
  integer :: my_ninterface
  integer :: nedges_coupled_loc
  integer :: nedges_acporo_coupled_loc
  integer :: nedges_elporo_coupled_loc

  character(len=256) :: prname


  do iproc = 0, nproc-1

    ! opens Database file
    write(prname, "('./OUTPUT_FILES/Database',i5.5)") iproc
    open(unit=15,file=trim(prname),status='unknown',iostat=ios)
    if( ios /= 0 ) stop 'error saving databases; check that directory OUTPUT_FILES exists'

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

    call write_glob2loc_nodes_database(15, iproc, npgeo, 1)

!   DK DK add support for using pml in mpi mode with external mesh
!   call write_partition_database(15, iproc, nspec, num_material, ngnod, 1)
    call write_partition_database(15, iproc, nspec, num_material, region_pml_external_mesh, ngnod, 1)

    write(15,*) 'npgeo nproc'
    write(15,*) npgeo,nproc

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

    write(15,*) 'initialfield add_Bielak_conditions'
    write(15,*) initialfield,add_Bielak_conditions

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

    write(15,*) 'save_ASCII_kernels'
    write(15,*) save_ASCII_kernels

    write(15,*) 'DRAW_SOURCES_AND_RECEIVERS'
    write(15,*) DRAW_SOURCES_AND_RECEIVERS

    write(15,*) 'Q0 freq0'
    write(15,*) Q0,freq0

    write(15,*) 'p_sv'
    write(15,*) p_sv

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

    write(15,*) 'ADD_SPRING_TO_STACEY'
    write(15,*) ADD_SPRING_TO_STACEY

    write(15,*) 'ADD_PERIODIC_CONDITIONS'
    write(15,*) ADD_PERIODIC_CONDITIONS

    write(15,*) 'PERIODIC_HORIZ_DIST'
    write(15,*) PERIODIC_HORIZ_DIST

    write(15,*) 'GPU_MODE'
    write(15,*) GPU_MODE

    write(15,*) 'nt deltat'
    write(15,*) nt,deltat

    write(15,*) 'NT_DUMP_ATTENUATION'
    write(15,*) NT_DUMP_ATTENUATION

    write(15,*) 'ACOUSTIC_FORCING'
    write(15,*) ACOUSTIC_FORCING

    write(15,*) 'NSOURCES'
    write(15,*) NSOURCES

    do i_source=1,NSOURCES
      write(15,*) 'source', i_source
      write(15,*) source_type(i_source),time_function_type(i_source)
      write(15,'(a100)') name_of_source_file(i_source) ! aXXX: Write wrong character if XXX != 100 !!!
      write(15,*) burst_band_width(i_source), &
                  xs(i_source),zs(i_source),f0(i_source),tshift_src(i_source), &
                  factor(i_source),anglesource(i_source), &
                  Mxx(i_source),Mzz(i_source),Mxz(i_source)
    enddo

    write(15,*) 'attenuation'
    write(15,*) N_SLS, f0_attenuation, READ_VELOCITIES_AT_f0

    write(15,*) 'Coordinates of macrobloc mesh (coorg):'

    call write_glob2loc_nodes_database(15, iproc, npgeo, 2)

    write(15,*) 'numat ngnod nspec pointsdisp plot_lowerleft_corner_only'
    write(15,*) nb_materials,ngnod,nspec,pointsdisp,plot_lowerleft_corner_only

    if (any_abs) then
      call write_abs_merge_database(15, iproc, 1)
    else
      nelemabs_loc = 0
    endif

    if (ACOUSTIC_FORCING) then
      call write_acoustic_forcing_merge_database(15, iproc, 1)
    else
      nelemacforcing_loc = 0
    endif

    call write_surface_database(15, nelem_acoustic_surface, acoustic_surface, nelem_acoustic_surface_loc, &
                              iproc, 1)

    call write_fluidsolid_edges_database(15,nedges_coupled, nedges_coupled_loc, &
                                        edges_coupled, iproc, 1)
    call write_fluidsolid_edges_database(15, nedges_acporo_coupled, nedges_acporo_coupled_loc, &
                                        edges_acporo_coupled, iproc, 1)
    call write_fluidsolid_edges_database(15, nedges_elporo_coupled, nedges_elporo_coupled_loc, &
                                        edges_elporo_coupled, iproc, 1)

    call write_axial_elements_database(15, nelem_on_the_axis, ispec_of_axial_elements,nelem_on_the_axis_loc,iproc,1, &
                                        remove_min_to_start_at_zero)

    if (.not. ( force_normal_to_surface .or. rec_normal_to_surface ) ) then
      nnodes_tangential_curve = 0
    endif

    write(15,*) 'nelemabs nelemacforcing nelem_acoustic_surface num_fluid_solid_edges'
    write(15,*) 'num_fluid_poro_edges num_solid_poro_edges'
    write(15,*) 'nnodes_tangential_curve nelem_on_the_axis'
    write(15,*) nelemabs_loc,nelemacforcing_loc,nelem_acoustic_surface_loc, &
                nedges_coupled_loc,nedges_acporo_coupled_loc,&
                nedges_elporo_coupled_loc,nnodes_tangential_curve, &
                nelem_on_the_axis_loc

    write(15,*) 'Material sets (num 1 rho vp vs 0 0 QKappa Qmu 0 0 0 0 0 0) or '
    write(15,*) '(num 2 rho c11 c13 c15 c33 c35 c55 c12 c23 c25 0 0 0) or '
    write(15,*) '(num 3 rhos rhof phi c k_xx k_xz k_zz Ks Kf Kfr etaf mufr Qmu)'
    do i=1,nb_materials
      if (icodemat(i) == ISOTROPIC_MATERIAL) then
         write(15,*) i,icodemat(i),rho_s(i),cp(i),cs(i),0,0,QKappa(i),Qmu(i),0,0,0,0,0,0
      else if (icodemat(i) == ANISOTROPIC_MATERIAL) then
         write(15,*) i,icodemat(i),rho_s(i), &
                    aniso3(i),aniso4(i),aniso5(i),aniso6(i),&
                    aniso7(i),aniso8(i),aniso9(i),aniso10(i),aniso11(i),aniso12(i),0,0
      else if(icodemat(i) == POROELASTIC_MATERIAL) then
         write(15,*) i,icodemat(i),rho_s(i),rho_f(i),phi(i),tortuosity(i), &
                    permxx(i),permxz(i),permzz(i),kappa_s(i),&
                    kappa_f(i),kappa_fr(i),eta_f(i),mu_fr(i),Qmu(i)
      else if (icodemat(i) <= 0) then ! The values will be read from an external tomo file
         write(15,*) i,icodemat(i),rho_s(i),cp(i),cs(i),0,0,QKappa(i),Qmu(i),0,0,0,0,0,0
      else
        stop 'Unknown material code'
      endif
    enddo

    write(15,*) 'Arrays kmato and knods for each bloc:'

!   DK DK add support for using PML in MPI mode with external mesh
!   call write_partition_database(15, iproc, nspec, num_material, ngnod, 2)
    call write_partition_database(15, iproc, nspec, num_material, region_pml_external_mesh, ngnod, 2)

    if ( nproc /= 1 ) then
      call write_interfaces_database(15, nproc, iproc, &
                              my_ninterface, my_interfaces, my_nb_interfaces, 1)

      write(15,*) 'Interfaces:'
      write(15,*) my_ninterface, maxval(my_nb_interfaces)

      call write_interfaces_database(15, nproc, iproc, &
                              my_ninterface, my_interfaces, my_nb_interfaces, 2)

    else
      write(15,*) 'Interfaces:'
      write(15,*) 0, 0
    endif


    write(15,*) 'List of absorbing elements (edge1 edge2 edge3 edge4 type):'
    if ( any_abs ) then
      call write_abs_merge_database(15, iproc, 2)
    endif

    write(15,*) 'List of acoustic forcing elements (edge1 edge2 edge3 edge4 type):'
    if ( ACOUSTIC_FORCING ) then
      call write_acoustic_forcing_merge_database(15, iproc, 2)
    endif

    write(15,*) 'List of acoustic free-surface elements:'
    call write_surface_database(15, nelem_acoustic_surface, acoustic_surface, nelem_acoustic_surface_loc, &
                                iproc, 2)


    write(15,*) 'List of acoustic elastic coupled edges:'
    call write_fluidsolid_edges_database(15, nedges_coupled, nedges_coupled_loc, &
                                        edges_coupled, iproc, 2)

    write(15,*) 'List of acoustic poroelastic coupled edges:'
    call write_fluidsolid_edges_database(15, nedges_acporo_coupled, nedges_acporo_coupled_loc, &
                                        edges_acporo_coupled, iproc, 2)

    write(15,*) 'List of poroelastic elastic coupled edges:'
    call write_fluidsolid_edges_database(15, nedges_elporo_coupled, nedges_elporo_coupled_loc, &
                                        edges_elporo_coupled, iproc, 2)

    write(15,*) 'List of tangential detection curve nodes:'
    !write(15,*) nnodes_tangential_curve
    write(15,*) force_normal_to_surface,rec_normal_to_surface

    if (force_normal_to_surface .or. rec_normal_to_surface) then
      do i = 1, nnodes_tangential_curve
        write(15,*) nodes_tangential_curve(1,i),nodes_tangential_curve(2,i)
      enddo
    endif

    write(15,*) 'List of axial elements:'
    call write_axial_elements_database(15, nelem_on_the_axis, ispec_of_axial_elements, nelem_on_the_axis_loc, iproc, 2, &
                                       remove_min_to_start_at_zero)

    ! closes Database file
    close(15)

  enddo

  end subroutine save_databases

