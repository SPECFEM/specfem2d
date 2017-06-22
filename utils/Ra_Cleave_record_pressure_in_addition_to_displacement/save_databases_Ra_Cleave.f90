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


  subroutine save_databases(nspec,num_material,region_pml_external_mesh, &
                            my_interfaces,my_nb_interfaces, &
                            nnodes_tangential_curve,nodes_tangential_curve )


! generates the databases for the solver

  use part_unstruct
  use parameter_file
  use source_file
  implicit none
  include "constants.h"

  integer :: nspec
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
    if ( ios /= 0 ) stop 'error saving databases; check that directory OUTPUT_FILES exists'

    write(15,*) '#'
    write(15,*) '# Database for SPECFEM2D'
    write(15,*) '# Dimitri Komatitsch, (c) University of Pau, France'
    write(15,*) '#'

    write(15,*) 'Title of the simulation'
    write(15,"(a100)") title

    write(15,*) 'Type of simulation'
    write(15,*) SIMULATION_TYPE, NOISE_TOMOGRAPHY, SAVE_FORWARD

    call write_glob2loc_nodes_database(15, iproc, npgeo, 1)

!   DK DK add support for using PML in MPI mode with external mesh
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

    write(15,*) 'seismotype seismo_p imagetype_postscript'
    write(15,*) seismotype,seismo_p,imagetype_postscript

    write(15,*) 'assign_external_model READ_EXTERNAL_SEP_FILE'
    write(15,*) assign_external_model,READ_EXTERNAL_SEP_FILE

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

    write(15,*) 'DRAW_SOURCES_AND_RECEIVERS'
    write(15,*) DRAW_SOURCES_AND_RECEIVERS

    write(15,*) 'Q0 freq0'
    write(15,*) Q0,freq0

    write(15,*) 'P_SV'
    write(15,*) P_SV

    write(15,*) 'factor_subsample_image'
    write(15,*) factor_subsample_image

    write(15,*) 'USE_SNAPSHOT_NUMBER_IN_FILENAME'
    write(15,*) USE_SNAPSHOT_NUMBER_IN_FILENAME

    write(15,*) 'DRAW_WATER_IN_BLUE'
    write(15,*) DRAW_WATER_IN_BLUE

    write(15,*) 'US_LETTER'
    write(15,*) US_LETTER

    write(15,*) 'POWER_DISPLAY_COLOR'
    write(15,*) POWER_DISPLAY_COLOR

    write(15,*) 'PERFORM_CUTHILL_MCKEE'
    write(15,*) PERFORM_CUTHILL_MCKEE

    write(15,*) 'SU_FORMAT'
    write(15,*) SU_FORMAT

    write(15,*) 'USER_T0'
    write(15,*) USER_T0

    write(15,*) 'time_stepping_scheme'
    write(15,*) time_stepping_scheme

    write(15,*) 'ADD_PERIODIC_CONDITIONS'
    write(15,*) ADD_PERIODIC_CONDITIONS

    write(15,*) 'PERIODIC_horiz_dist'
    write(15,*) PERIODIC_horiz_dist

    write(15,*) 'PERIODIC_DETECT_TOL'
    write(15,*) PERIODIC_DETECT_TOL

    write(15,*) 'nt deltat'
    write(15,*) nt,deltat
    write(15,*) 'NSOURCES'
    write(15,*) NSOURCES

    do i_source=1,NSOURCES
      write(15,*) 'source', i_source
      write(15,*) source_type(i_source),time_function_type(i_source), &
                  xs(i_source),zs(i_source),f0_source(i_source),tshift_src(i_source), &
                  factor(i_source),anglesource(i_source), &
                  Mxx(i_source),Mzz(i_source),Mxz(i_source)
    enddo

    write(15,*) 'attenuation'
    write(15,*) N_SLS, f0_attenuation

    write(15,*) 'Coordinates of macrobloc mesh (coorg):'

    call write_glob2loc_nodes_database(15, iproc, npgeo, 2)

    write(15,*) 'numat ngnod nspec pointsdisp plot_lowerleft_corner_only'
    write(15,*) nb_materials,ngnod,nspec,pointsdisp,plot_lowerleft_corner_only

    if (any_abs) then
      call write_abs_merge_database(15, iproc, 1)
    else
      nelemabs_loc = 0
    endif

    call write_surface_database(15, nelem_acoustic_surface, acoustic_surface, nelem_acoustic_surface_loc, &
                              iproc, 1)

    call write_fluidsolid_edges_database(15,nedges_coupled, nedges_coupled_loc, &
                                        edges_coupled, iproc, 1)
    call write_fluidsolid_edges_database(15, nedges_acporo_coupled, nedges_acporo_coupled_loc, &
                                        edges_acporo_coupled, iproc, 1)
    call write_fluidsolid_edges_database(15, nedges_elporo_coupled, nedges_elporo_coupled_loc, &
                                        edges_elporo_coupled, iproc, 1)

    if (.not. ( force_normal_to_surface .or. rec_normal_to_surface ) ) then
      nnodes_tangential_curve = 0
    endif

    write(15,*) 'nelemabs nelem_acoustic_surface num_fluid_solid_edges num_fluid_poro_edges'
    write(15,*) 'num_solid_poro_edges nnodes_tangential_curve'
    write(15,*) nelemabs_loc,nelem_acoustic_surface_loc, &
                nedges_coupled_loc,nedges_acporo_coupled_loc, &
                nedges_elporo_coupled_loc,nnodes_tangential_curve

    write(15,*) 'Material sets (num 1 rho vp vs 0 0 QKappa Qmu 0 0 0 0 0 0) or '
    write(15,*) '(num 2 rho c11 c13 c33 c44 QKappa Qmu 0 0 0 0 0 0) or '
    write(15,*) '(num 3 rhos rhof phi c k_xx k_xz k_zz Ks Kf Kfr etaf mufr Qmu)'
    do i=1,nb_materials
      if (icodemat(i) == ISOTROPIC_MATERIAL) then
         write(15,*) i,icodemat(i),rho_s(i),cp(i),cs(i),0,0,QKappa(i),Qmu(i),0,0,0,0,0,0
      else if (icodemat(i) == POROELASTIC_MATERIAL) then
         write(15,*) i,icodemat(i),rho_s(i),rho_f(i),phi(i),tortuosity(i), &
                    permxx(i),permxz(i),permzz(i),kappa_s(i), &
                    kappa_f(i),kappa_fr(i),eta_f(i),mu_fr(i),Qmu(i)
      else
         write(15,*) i,icodemat(i),rho_s(i),cp(i),cs(i), &
                    aniso3(i),aniso4(i),aniso5(i),aniso6(i), &
                    aniso7(i),aniso8(i),QKappa(i),Qmu(i),0,0
      endif
    enddo

    write(15,*) 'Arrays kmato and knods for each bloc:'

!   DK DK add support for using PML in MPI mode with external mesh
!   call write_partition_database(15, iproc, nspec, num_material, ngnod, 2)
    call write_partition_database(15, iproc, nspec, num_material, region_pml_external_mesh, ngnod, 2)

    if (nproc /= 1) then
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
    if (any_abs) then
      call write_abs_merge_database(15, iproc, 2)
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

    ! closes Database file
    close(15)

  enddo

  end subroutine save_databases

