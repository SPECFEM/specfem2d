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

  subroutine compute_forces_viscoacoustic_GPU(compute_b_wavefield)

  use specfem_par, only: NPROC,ninterface,max_nibool_interfaces_ext_mesh,nibool_interfaces_ext_mesh, &
                         my_neighbors,ninterface_acoustic,inum_interfaces_acoustic, &
                         nelem_acoustic_surface,num_fluid_solid_edges,UNDO_ATTENUATION_AND_OR_PML, &
                         STACEY_ABSORBING_CONDITIONS,PML_BOUNDARY_CONDITIONS, &
                         coupled_acoustic_elastic,coupled_acoustic_poro, &
                         SIMULATION_TYPE,ATTENUATION_VISCOACOUSTIC, &
                         deltatover2,b_deltatover2

  use specfem_par, only: nspec_outer_acoustic, nspec_inner_acoustic,NO_BACKWARD_RECONSTRUCTION

  use specfem_par_gpu, only: Mesh_pointer, &
                             buffer_send_scalar_gpu,buffer_recv_scalar_gpu, &
                             b_buffer_send_scalar_gpu,b_buffer_recv_scalar_gpu, &
                             request_send_recv_scalar_gpu,b_request_send_recv_scalar_gpu

  implicit none

  ! parameter useful only for undo_attenuation (to know which wavefield to compute)
  logical,intent(in) :: compute_b_wavefield

  ! local parameters
  integer :: iphase
  logical :: compute_wavefield_1    ! forward wavefield (forward or adjoint for SIMULATION_TYPE > 1)
  logical :: compute_wavefield_2    ! backward/reconstructed wavefield (b_** arrays)

  ! determines which wavefields to compute
  if ((.not. UNDO_ATTENUATION_AND_OR_PML) .and. (SIMULATION_TYPE == 1 .or. NO_BACKWARD_RECONSTRUCTION) ) then
    ! forward wavefield only
    compute_wavefield_1 = .true.
    compute_wavefield_2 = .false.
  else if ((.not. UNDO_ATTENUATION_AND_OR_PML) .and. SIMULATION_TYPE == 3) then
    ! forward & backward wavefields
    compute_wavefield_1 = .true.
    compute_wavefield_2 = .true.
  else if (UNDO_ATTENUATION_AND_OR_PML .and. compute_b_wavefield) then
    ! only backward wavefield
    compute_wavefield_1 = .false.
    compute_wavefield_2 = .true.
  else
    ! default forward wavefield only
    compute_wavefield_1 = .true.
    compute_wavefield_2 = .false.
  endif

  ! coupled simulation
  ! requires different coupling terms for forward/adjoint and backpropagated wavefields
  if (coupled_acoustic_elastic) then
    if (SIMULATION_TYPE == 3) then
      if (compute_b_wavefield) then
        ! only backward wavefield
        compute_wavefield_1 = .false.
        compute_wavefield_2 = .true.
      else
        ! only forward/adjoint wavefield
        compute_wavefield_1 = .true.
        compute_wavefield_2 = .false.
      endif
    endif
  endif

  if (nelem_acoustic_surface > 0) then
    !  enforces free surface (zeroes potentials at free surface)
    call acoustic_enforce_free_surf_cuda(Mesh_pointer,compute_wavefield_1,compute_wavefield_2)
  endif

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase = 1,2
    ! acoustic pressure term
    ! includes code for SIMULATION_TYPE==3
    call compute_forces_acoustic_cuda(Mesh_pointer, iphase, &
                                      nspec_outer_acoustic, nspec_inner_acoustic,ATTENUATION_VISCOACOUSTIC, &
                                      compute_wavefield_1,compute_wavefield_2)

    ! PML boundary conditions enforces zero potentials on boundary
    if (PML_BOUNDARY_CONDITIONS) then
      call pml_boundary_acoustic_cuda(Mesh_pointer,compute_wavefield_1,compute_wavefield_2)
    endif

    ! computes additional contributions
    if (iphase == 1) then
      ! Stacey absorbing boundary conditions
      if (STACEY_ABSORBING_CONDITIONS) then
        call compute_stacey_acoustic_GPU(iphase,compute_wavefield_1,compute_wavefield_2)
      endif

      ! elastic coupling
      if (coupled_acoustic_elastic) then
        if (compute_wavefield_1) &
          call compute_coupling_ac_el_cuda(Mesh_pointer,iphase,num_fluid_solid_edges,1) ! 1 == forward/adjoint
        if (compute_wavefield_2) &
          call compute_coupling_ac_el_cuda(Mesh_pointer,iphase,num_fluid_solid_edges,3) ! 3 == backward
      endif

      ! poroelastic coupling
      if (coupled_acoustic_poro) then
        call stop_the_code('poroelastic not implemented yet in GPU mode')
      endif

      ! sources
      call compute_add_sources_acoustic_GPU(iphase,compute_wavefield_1,compute_wavefield_2)
    endif

    ! assemble all the contributions between slices using MPI
    if (ninterface_acoustic > 0) then

      if (iphase == 1) then
        if (compute_wavefield_1) then
          ! sends potential_dot_dot_acoustic values to corresponding MPI interface neighbors (non-blocking)
          call transfer_boun_pot_from_device(Mesh_pointer, &
                                             buffer_send_scalar_gpu, &
                                             1) ! -- 1 == fwd accel

          call assemble_MPI_scalar_send_cuda(NPROC, &
                                             buffer_send_scalar_gpu,buffer_recv_scalar_gpu, &
                                             ninterface,max_nibool_interfaces_ext_mesh, &
                                             nibool_interfaces_ext_mesh, &
                                             my_neighbors, &
                                             request_send_recv_scalar_gpu,ninterface_acoustic,inum_interfaces_acoustic)
        endif
        ! adjoint simulations
        if (compute_wavefield_2) then
          call transfer_boun_pot_from_device(Mesh_pointer, &
                                             b_buffer_send_scalar_gpu, &
                                             3) ! -- 3 == adjoint b_accel

          call assemble_MPI_scalar_send_cuda(NPROC, &
                                             b_buffer_send_scalar_gpu,b_buffer_recv_scalar_gpu, &
                                             ninterface,max_nibool_interfaces_ext_mesh, &
                                             nibool_interfaces_ext_mesh, &
                                             my_neighbors, &
                                             b_request_send_recv_scalar_gpu,ninterface_acoustic,inum_interfaces_acoustic)
        endif

      else
        if ( compute_wavefield_1) &
          ! waits for send/receive requests to be completed and assembles values
          call assemble_MPI_scalar_write_cuda(NPROC, &
                                              Mesh_pointer, &
                                              buffer_recv_scalar_gpu, &
                                              ninterface, &
                                              max_nibool_interfaces_ext_mesh, &
                                              request_send_recv_scalar_gpu, &
                                              1,ninterface_acoustic,inum_interfaces_acoustic)



        ! adjoint simulations
        if (compute_wavefield_2) then
          call assemble_MPI_scalar_write_cuda(NPROC, &
                                              Mesh_pointer, &
                                              b_buffer_recv_scalar_gpu, &
                                              ninterface, &
                                              max_nibool_interfaces_ext_mesh, &
                                              b_request_send_recv_scalar_gpu, &
                                              3,ninterface_acoustic,inum_interfaces_acoustic)
        endif
      endif !iphase

    endif !interface_acoustic

  enddo !iphase

  ! TODO PML saves interface values, record the field at each time step on the PML edges and write it to file.
  ! This would have to be done here (in GPU mode) also if one wants to do adjoint simulations
  ! without using NO_BACKWARD_RECONSTRUCTION
  !if (PML_BOUNDARY_CONDITIONS) then
  !  if (nglob_interface > 0) then
  !    if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
  !      do i = 1, nglob_interface
  !        write(72) potential_dot_dot_acoustic(point_interface(i)), potential_dot_acoustic(point_interface(i)), &
  !                  potential_acoustic(point_interface(i))
  !      enddo
  !    endif
  !  endif
  !endif

  ! divides pressure with mass matrix
  ! corrector:
  ! updates the chi_dot term which requires chi_dot_dot(t+delta)
  call kernel_3_acoustic_cuda(Mesh_pointer,deltatover2,b_deltatover2,compute_wavefield_1,compute_wavefield_2)

  ! enforces free surface (zeroes potentials at free surface)
  call acoustic_enforce_free_surf_cuda(Mesh_pointer,compute_wavefield_1,compute_wavefield_2)

  end subroutine compute_forces_viscoacoustic_GPU

!
!---------------------------------------------------------------------------------------------
!

! for acoustic solver on GPU

  subroutine compute_stacey_acoustic_GPU(iphase,compute_wavefield_1,compute_wavefield_2)

  use constants, only: CUSTOM_REAL,NGLLX

  use specfem_par, only: nspec_bottom,nspec_left,nspec_top,nspec_right, &
                         b_absorb_acoustic_left,b_absorb_acoustic_right, &
                         b_absorb_acoustic_bottom, b_absorb_acoustic_top, &
                         it,NSTEP,SIMULATION_TYPE,SAVE_FORWARD, &
                         num_abs_boundary_faces,UNDO_ATTENUATION_AND_OR_PML,NO_BACKWARD_RECONSTRUCTION

  use specfem_par_gpu, only: Mesh_pointer

  implicit none

  ! communication overlap
  integer,intent(in) :: iphase
  logical,intent(in) :: compute_wavefield_1,compute_wavefield_2

  ! adjoint simulations
  real(kind=CUSTOM_REAL),dimension(NGLLX,nspec_bottom) :: b_absorb_potential_bottom_slice
  real(kind=CUSTOM_REAL),dimension(NGLLX,nspec_left) :: b_absorb_potential_left_slice
  real(kind=CUSTOM_REAL),dimension(NGLLX,nspec_right) :: b_absorb_potential_right_slice
  real(kind=CUSTOM_REAL),dimension(NGLLX,nspec_top) :: b_absorb_potential_top_slice

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return
  if (NO_BACKWARD_RECONSTRUCTION) return

  if (SIMULATION_TYPE == 3 .and. (.not. UNDO_ATTENUATION_AND_OR_PML) ) then
    ! gets absorbing contributions buffers
    b_absorb_potential_left_slice(:,:) = b_absorb_acoustic_left(:,:,NSTEP-it+1)
    b_absorb_potential_right_slice(:,:) = b_absorb_acoustic_right(:,:,NSTEP-it+1)
    b_absorb_potential_top_slice(:,:) = b_absorb_acoustic_top(:,:,NSTEP-it+1)
    b_absorb_potential_bottom_slice(:,:) = b_absorb_acoustic_bottom(:,:,NSTEP-it+1)
  endif

  ! absorbs absorbing-boundary surface using Sommerfeld condition (vanishing field in the outer-space)
  call compute_stacey_acoustic_cuda(Mesh_pointer,iphase, &
                                    b_absorb_potential_left_slice,b_absorb_potential_right_slice, &
                                    b_absorb_potential_top_slice,b_absorb_potential_bottom_slice, &
                                    compute_wavefield_1,compute_wavefield_2,UNDO_ATTENUATION_AND_OR_PML)

  ! adjoint simulations: stores absorbed wavefield part
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. (.not. UNDO_ATTENUATION_AND_OR_PML)) then
    ! writes out absorbing boundary value only when second phase is running
    b_absorb_acoustic_bottom(:,:,it) = b_absorb_potential_bottom_slice(:,:)
    b_absorb_acoustic_right(:,:,it) = b_absorb_potential_right_slice(:,:)
    b_absorb_acoustic_top(:,:,it) = b_absorb_potential_top_slice(:,:)
    b_absorb_acoustic_left(:,:,it) = b_absorb_potential_left_slice(:,:)
  endif

  end subroutine compute_stacey_acoustic_GPU

!
!---------------------------------------------------------------------------------------------
!

  subroutine compute_add_sources_acoustic_GPU(iphase,compute_wavefield_1,compute_wavefield_2)

  use specfem_par, only: it,SIMULATION_TYPE,NSTEP,nadj_rec_local, &
    SOURCE_IS_MOVING,NSOURCES

  use moving_sources_par, only: nsources_local_moving

  use specfem_par_gpu, only: Mesh_pointer

  implicit none

  ! communication overlap
  integer,intent(in) :: iphase
  logical,intent(in) :: compute_wavefield_1,compute_wavefield_2

  ! local parameters
  integer :: it_tmp

  ! forward simulations
  if (SIMULATION_TYPE == 1) then
    if (SOURCE_IS_MOVING) then
      ! sources move position
      call compute_add_moving_sources_ac_cuda(Mesh_pointer, iphase, nsources_local_moving(it), it, NSTEP, NSOURCES)

      ! Old version. Using this version you can comment call init_moving_sources() in prepare_gpu.f90
      ! call compute_add_sources_acoustic_GPU_moving_sources_old(iphase, compute_wavefield_1, compute_wavefield_2)
    else
      ! sources stay at same location
      call compute_add_sources_ac_cuda(Mesh_pointer,iphase,it)
    endif
  endif

  ! adjoint simulations
  ! time step index
  it_tmp = NSTEP - it + 1

  ! adds acoustic adjoint sources
  if (SIMULATION_TYPE /= 1 .and. nadj_rec_local > 0 .and. compute_wavefield_1) &
    call add_sources_ac_sim_2_or_3_cuda(Mesh_pointer, iphase, it_tmp, nadj_rec_local, NSTEP)

  ! kernel simulations w/ backward/reconstructed sources
  if (compute_wavefield_2) call compute_add_sources_ac_s3_cuda(Mesh_pointer, iphase, it_tmp)

  end subroutine compute_add_sources_acoustic_GPU

!
!---------------------------------------------------------------------------------------------
!

!unused routine, left here for reference
!
!    subroutine compute_add_sources_acoustic_GPU_moving_sources_old(iphase,compute_wavefield_1,compute_wavefield_2)
!    ! This subroutine is useful to get an idea of how moving sources are handled
!    ! ... but it is very slow !! (a lot of communication with the GPU)
!    ! See new version in prepare_gpu.f90
!
!    use constants, only: NGLLX,NGLLZ,TINYVAL
!
!    use specfem_par
!    use specfem_par_gpu
!
!    implicit none
!
!    ! communication overlap
!    integer,intent(in) :: iphase
!    logical,intent(in) :: compute_wavefield_1,compute_wavefield_2
!
!    !local variables
!    integer :: i_source,i,j,k,ispec
!    double precision :: hlagrange
!    double precision :: xInitSource,zInitSource,xsrc,zsrc,vxSource,vzSource, &
!                        time_val,t_used
!    ! single source array
!    real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: sourcearray
!
!    ! forward simulations
!    if (SIMULATION_TYPE == 1) then
!      xInitSource = x_source(i_source) ! m
!      zInitSource = z_source(i_source) ! m
!      vxSource = vx_source(i_source) ! m/s
!      vzSource = vz_source(i_source) ! m/s
!
!      if (time_stepping_scheme == 1) then
!        ! Newmark
!        time_val = (it-1)*deltat
!      else
!        call exit_MPI(myrank,'Only Newmark time scheme is implemented for moving sources (0)')
!      endif
!
!      ! moves and re-locates sources along x-axis
!      do i_source = 1,NSOURCES
!        if (abs(source_time_function(i_source,it,i_stage)) > TINYVAL) then
!          t_used = (time_val-t0-tshift_src(i_source))
!
!          xsrc = xInitSource + vxSource*t_used
!          zsrc = zInitSource + vzSource*t_used
!
!          ! collocated force source
!          call locate_source(ibool,coord,nspec,nglob,xigll,zigll, &
!                             xsrc,zsrc, &
!                             ispec_selected_source(i_source),islice_selected_source(i_source), &
!                             NPROC,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,NGNOD,npgeo, &
!                             iglob_source(i_source),.true.)
!
!          call lagrange_any(xi_source(i_source),NGLLX,xigll,hxis,hpxis)
!          call lagrange_any(gamma_source(i_source),NGLLZ,zigll,hgammas,hpgammas)
!
!          ! stores Lagrangians for source
!          hxis_store(i_source,:) = hxis(:)
!          hgammas_store(i_source,:) = hgammas(:)
!
!          sourcearray(:,:,:) = 0._CUSTOM_REAL
!
!          ! element containing source
!          ispec = ispec_selected_source(i_source)
!
!          if (myrank == islice_selected_source(i_source)) then
!
!            ! computes source arrays
!            if (source_type(i_source) == 1) then
!              ! collocated force source
!              do j = 1,NGLLZ
!                do i = 1,NGLLX
!                  hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
!                  ! source element is acoustic
!                  if (ispec_is_acoustic(ispec)) then
!                    ! sourcearray(1 and 2, i, j) are set to hlagrange:
!                    sourcearray(:,i,j) = hlagrange
!                  else
!                    call exit_MPI(myrank,'Moving source not in acoustic element (GPU, not implemented)')
!                  endif
!                enddo
!              enddo
!            else
!              call exit_MPI(myrank,'Moving source with source_type != 1 (GPU, not implemented)')
!            endif
!
!            ! stores sourcearray for all sources
!            sourcearrays(:,:,:,i_source) = sourcearray(:,:,:)
!
!          endif
!        endif
!      enddo
!
!      ! counts sources in this process slice
!      nsources_local = 0
!      do i = 1, NSOURCES
!        if (myrank == islice_selected_source(i)) then
!          nsources_local = nsources_local + 1
!        endif
!      enddo
!      deallocate(ispec_selected_source_loc)
!      allocate(ispec_selected_source_loc(nsources_local))
!      ispec_selected_source_loc(:) = 1
!      j = 0
!      do i = 1, NSOURCES
!        if (myrank == islice_selected_source(i)) then
!          if (j > nsources_local) call stop_the_code('Error with the number of local sources')
!          j = j + 1
!          ispec_selected_source_loc(j)  = ispec_selected_source(i)
!        endif
!      enddo
!      deallocate(sourcearray_loc)
!      if (nsources_local > 0) then
!        allocate(sourcearray_loc(nsources_local,NDIM,NGLLX,NGLLX))
!      else
!        allocate(sourcearray_loc(1,1,1,1))
!      endif
!
!      sourcearray_loc(:,:,:,:) = 0._CUSTOM_REAL
!      k = 0
!      do i_source = 1,NSOURCES
!        if (myrank == islice_selected_source(i_source)) then
!          ! source belongs to this process
!          k = k + 1
!          sourcearray_loc(k,:,:,:) = sourcearrays(:,:,:,i_source)
!        endif
!      enddo
!
!      ! Send these values to device:
!      call recompute_source_position_cuda(Mesh_pointer, nsources_local, sourcearray_loc, &
!                                      ispec_selected_source_loc)
!
!      call compute_add_sources_ac_cuda(Mesh_pointer, iphase, it)
!    endif
!    ! adjoint simulations
!    if (SIMULATION_TYPE == 3 .and. nadj_rec_local > 0 .and. compute_wavefield_1) &
!      call add_sources_ac_sim_2_or_3_cuda(Mesh_pointer,iphase, NSTEP -it + 1, nadj_rec_local,NSTEP)
!
!    ! adjoint simulations
!    if (compute_wavefield_2) call compute_add_sources_ac_s3_cuda(Mesh_pointer, iphase, NSTEP -it + 1)
!
!    end subroutine compute_add_sources_acoustic_GPU_moving_sources_old
!
