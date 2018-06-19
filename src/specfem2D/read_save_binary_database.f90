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

  subroutine save_binary_database()

  use constants
  use specfem_par
  use specfem_par_gpu

  implicit none

  integer ier
  character(len=MAX_STRING_LEN) :: outputname


  ! saves data in a binary file
  write(outputname,'(a,i6.6,a)') 'proc',myrank,'_data.bin'
        open(unit = 2040, file = trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
        if (ier /= 0) call stop_the_code('Error writing data file to disk')

  write(2040) nglob,ELASTIC_SIMULATION,POROELASTIC_SIMULATION, &
              ACOUSTIC_SIMULATION,coupled_acoustic_elastic, &
              any_acoustic,any_elastic,any_poroelastic

  write(2040) coord,jacobian,xix,xiz,gammax,gammaz, &
              this_ibool_is_a_periodic_edge,ibool,ispec_is_inner

  if (GPU_MODE) then
    write (2040) abs_boundary_ij,abs_boundary_normal,abs_boundary_jacobian1Dw, &
                 free_ac_ispec,cote_abs,free_surface_ij,ANY_ANISOTROPY
  endif

  if (any_acoustic) then
    write(2040) rmass_inverse_acoustic,num_phase_ispec_acoustic
    write(2040) phase_ispec_inner_acoustic,nspec_inner_acoustic,nspec_outer_acoustic,ispec_is_acoustic
    if (ATTENUATION_VISCOACOUSTIC) write(2040) rmass_inverse_e1
  endif

  if (any_elastic) then
    write(2040) rmass_inverse_elastic,num_phase_ispec_elastic
    write(2040) phase_ispec_inner_elastic,nspec_inner_elastic,nspec_outer_elastic,ispec_is_elastic
  endif

  if (coupled_acoustic_elastic) then
    write(2040) fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge,fluid_solid_elastic_ispec,fluid_solid_elastic_iedge, &
                ivalue_inverse,jvalue_inverse,ivalue,jvalue
    if (GPU_MODE) write(2040) coupling_ac_el_ispec,coupling_ac_el_ij,coupling_ac_el_normal,coupling_ac_el_jacobian1Dw
  endif

  if (NPROC > 1) then
    write(2040) ninterface_acoustic,ninterface_elastic,inum_interfaces_acoustic,inum_interfaces_elastic, &
                nibool_interfaces_acoustic,nibool_interfaces_elastic,nibool_interfaces_ext_mesh, &
                ibool_interfaces_acoustic,ibool_interfaces_elastic,ibool_interfaces_ext_mesh
  endif

  close(2040)

  ! saves all data regarding sources in a binary file
  write(outputname,'(a,i6.6,a)') 'proc',myrank,'_sources_info.bin'
      open(unit = 2040, file = trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing sources info file to disk')

  write(2040) source_time_function,nsources_local,sourcearrays,islice_selected_source,ispec_selected_source

  close(2040)

  ! saves all data regarding receivers in a binary file
  write(outputname,'(a,i6.6,a)') 'proc',myrank,'_receivers_info.bin'
      open(unit = 2040, file = trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing receivers info data file to disk')

  write(2040) nrecloc,nrec
  write(2040) recloc,ispec_selected_rec_loc,cosrot_irec,sinrot_irec,xir_store_loc,gammar_store_loc,st_xval,st_zval, &
              station_name,network_name,islice_selected_rec

  close(2040)

  end subroutine save_binary_database

!
!-----------------------------------------------------------------------------------
!

  subroutine read_binary_database_part1()

  use constants
  use specfem_par
  use specfem_par_gpu
  use specfem_par_movie

  implicit none

  integer ier
  character(len=MAX_STRING_LEN) :: outputname

  ! read setup data from a binary file that need to be known in advance (for allocations purpose)
  write(outputname,'(a,i6.6,a)') trim(OUTPUT_FILES)//'proc',myrank,'_data.bin'
      open(unit=2040,file=outputname,status='old',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_data.bin')

  read(2040) nglob,ELASTIC_SIMULATION,POROELASTIC_SIMULATION, &
             ACOUSTIC_SIMULATION,coupled_acoustic_elastic, &
             any_acoustic,any_elastic,any_poroelastic

  if (any_poroelastic) &
      call exit_MPI(myrank,'currently cannot have database mode if poroelastic simulation')
  if (PML_BOUNDARY_CONDITIONS) &
      call exit_MPI(myrank,'currently cannot have database mode with PML')
  if (AXISYM) &
      call exit_MPI(myrank,'currently cannot have database mode with AXISYM')
  if (output_postscript_snapshot) &
      call exit_MPI(myrank,'currently cannot have database mode with postscript')
  if (GPU_MODE) &
      call exit_MPI(myrank,'currently cannot have database mode with GPU_MODE')

  end subroutine read_binary_database_part1
!
!-----------------------------------------------------------------------------------
!

  subroutine read_binary_database_part2()

  use constants
  use specfem_par
  use specfem_par_gpu
  use specfem_par_movie

  implicit none

  integer :: ier,n_sls_loc

  ! reads setup data from a binary file
  allocate(ispec_is_inner(nspec))
  allocate(this_ibool_is_a_periodic_edge(nglob),coord(NDIM,nglob),stat=ier)

  read(2040) coord,jacobian,xix,xiz,gammax,gammaz, &
             this_ibool_is_a_periodic_edge,ibool,ispec_is_inner

  if (GPU_MODE) then
    read (2040) abs_boundary_ij,abs_boundary_normal,abs_boundary_jacobian1Dw, &
                free_ac_ispec,cote_abs,free_surface_ij,ANY_ANISOTROPY
  endif

  if (any_acoustic) then
    read(2040) rmass_inverse_acoustic,num_phase_ispec_acoustic

    allocate( phase_ispec_inner_acoustic(num_phase_ispec_acoustic,2),stat=ier)
    if (ier /= 0 ) call stop_the_code('Error allocating array phase_ispec_inner_acoustic')

    read(2040) phase_ispec_inner_acoustic,nspec_inner_acoustic,nspec_outer_acoustic,ispec_is_acoustic
  else
    ! allocates dummy array
    num_phase_ispec_acoustic = 0
    allocate( phase_ispec_inner_acoustic(num_phase_ispec_acoustic,2),stat=ier)
    if (ier /= 0 ) call stop_the_code('Error allocating dummy array phase_ispec_inner_acoustic')
  endif

  if (any_elastic) then
    read(2040) rmass_inverse_elastic,num_phase_ispec_elastic

    allocate( phase_ispec_inner_elastic(num_phase_ispec_elastic,2),stat=ier)
    if (ier /= 0 ) call stop_the_code('Error allocating array phase_ispec_inner_elastic')

    read(2040) phase_ispec_inner_elastic,nspec_inner_elastic,nspec_outer_elastic,ispec_is_elastic
  else
    ! allocates dummy array
    num_phase_ispec_elastic = 0
    allocate( phase_ispec_inner_elastic(num_phase_ispec_elastic,2),stat=ier)
    if (ier /= 0 ) call stop_the_code('Error allocating dummy array phase_ispec_inner_elastic')
  endif

  if (coupled_acoustic_elastic) then
    read(2040) fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge,fluid_solid_elastic_ispec,fluid_solid_elastic_iedge, &
               ivalue_inverse,jvalue_inverse,ivalue,jvalue
    if (GPU_MODE) read(2040) coupling_ac_el_ispec,coupling_ac_el_ij,coupling_ac_el_normal,coupling_ac_el_jacobian1Dw
  endif

  if (NPROC > 1) then
    read(2040) ninterface_acoustic,ninterface_elastic,inum_interfaces_acoustic,inum_interfaces_elastic, &
               nibool_interfaces_acoustic,nibool_interfaces_elastic,nibool_interfaces_ext_mesh, &
               ibool_interfaces_acoustic,ibool_interfaces_elastic,ibool_interfaces_ext_mesh

    max_nibool_interfaces_ext_mesh = maxval(nibool_interfaces_ext_mesh(:))
    max_ibool_interfaces_size_ac = maxval(nibool_interfaces_acoustic(:))
    max_ibool_interfaces_size_el = NDIM*maxval(nibool_interfaces_elastic(:))

    ! allocations
    if (ACOUSTIC_SIMULATION) then
      n_sls_loc = 0
      if (ATTENUATION_VISCOACOUSTIC) n_sls_loc = N_SLS
      allocate(request_send_recv_acoustic(ninterface_acoustic*2),stat=ier)
      if (ier /= 0) call stop_the_code('error in allocation of array request_send_recv_acoustic')
      allocate(buffer_send_faces_vector_ac(max_ibool_interfaces_size_ac*(n_sls_loc+1),ninterface_acoustic),stat=ier)
      if (ier /= 0) call stop_the_code('error in allocation of array buffer_send_faces_vector_ac')
      allocate(buffer_recv_faces_vector_ac(max_ibool_interfaces_size_ac*(n_sls_loc+1),ninterface_acoustic),stat=ier)
      if (ier /= 0) call stop_the_code('error in allocation of array buffer_recv_faces_vector_ac')
    endif

    if (ELASTIC_SIMULATION) then
      allocate(request_send_recv_elastic(ninterface_elastic*2),stat=ier)
      if (ier /= 0) call stop_the_code('error in allocation of array request_send_recv_elastic')
      allocate(buffer_send_faces_vector_el(max_ibool_interfaces_size_el,ninterface_elastic),stat=ier)
      if (ier /= 0) call stop_the_code('error in allocation of array buffer_send_faces_vector_el')
      allocate(buffer_recv_faces_vector_el(max_ibool_interfaces_size_el,ninterface_elastic),stat=ier)
      if (ier /= 0) call stop_the_code('error in allocation of array buffer_recv_faces_vector_el')
    endif

  endif

  close(2040)

  allocate(vector_field_display(NDIM,nglob))
  vector_field_display(:,:) = 0.d0

  end subroutine read_binary_database_part2

!
!-----------------------------------------------------------------------------------
!

  subroutine read_sources_receivers()

  use constants
  use specfem_par
  use specfem_par_gpu
  use specfem_par_movie

  implicit none

  integer ier
  character(len=MAX_STRING_LEN) :: outputname

  ! reads all data regarding sources from a binary file
  write(outputname,'(a,i6.6,a)') 'proc',myrank,'_sources_info.bin'
      open(unit = 2040, file = trim(OUTPUT_FILES)//outputname,status='old',action='read',form='unformatted', iostat=ier)
      if (ier /= 0) call stop_the_code('Error reading sources info from disk')


  if (initialfield) then
    allocate(source_time_function(1,1,1))
  else
    allocate(source_time_function(NSOURCES,NSTEP,stage_time_scheme),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating array source_time_function')
  endif

  allocate(hxis_store(NSOURCES,NGLLX), &
           hgammas_store(NSOURCES,NGLLZ),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating source h**_store arrays')

  read(2040) source_time_function,nsources_local,sourcearrays,islice_selected_source,ispec_selected_source

  close(2040)

  ! reads all data from receivers a binary file
  write(outputname,'(a,i6.6,a)') 'proc',myrank,'_receivers_info.bin'
      open(unit = 2040, file = trim(OUTPUT_FILES)//outputname,status='old',action='read',form='unformatted', iostat=ier)
      if (ier /= 0) call stop_the_code('Error reading receivers info from disk')

  read(2040) nrecloc,nrec

  allocate(st_xval(nrec), &
           st_zval(nrec), &
           station_name(nrec), &
           network_name(nrec), &
           recloc(nrec), &
           islice_selected_rec(nrec), &
           x_final_receiver(nrec), &
           z_final_receiver(nrec),stat=ier)

  allocate(ispec_selected_rec_loc(nrecloc))

  ! allocate Lagrange interpolators for receivers
  allocate(xir_store_loc(nrecloc,NGLLX), &
           gammar_store_loc(nrecloc,NGLLZ),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating local receiver h**_store arrays')

  allocate(anglerec_irec(nrecloc), &
           cosrot_irec(nrecloc), &
           sinrot_irec(nrecloc), &
           rec_tangential_detection_curve(nrecloc),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating tangential arrays')

  read(2040) recloc,ispec_selected_rec_loc,cosrot_irec,sinrot_irec,xir_store_loc,gammar_store_loc,st_xval,st_zval, &
             station_name,network_name,islice_selected_rec

  close(2040)

  end subroutine read_sources_receivers
