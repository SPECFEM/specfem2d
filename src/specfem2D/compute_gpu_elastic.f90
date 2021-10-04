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

! elastic solver

  subroutine compute_forces_viscoelastic_GPU(compute_b_wavefield)

  use constants, only: CUSTOM_REAL,NGLLX,NDIM

  use specfem_par, only: myrank,NPROC,ninterface,max_nibool_interfaces_ext_mesh,nibool_interfaces_ext_mesh, &
    my_neighbors,ninterface_elastic,inum_interfaces_elastic,ibool_interfaces_ext_mesh, &
    num_fluid_solid_edges,UNDO_ATTENUATION_AND_OR_PML, &
    STACEY_ABSORBING_CONDITIONS,PML_BOUNDARY_CONDITIONS, &
    coupled_acoustic_elastic,coupled_elastic_poro, &
    SIMULATION_TYPE,ATTENUATION_VISCOELASTIC, &
    deltat,deltatover2,b_deltatover2

  use specfem_par, only: nspec_outer_elastic,nspec_inner_elastic,any_anisotropy,NO_BACKWARD_RECONSTRUCTION

  use specfem_par_gpu, only: Mesh_pointer, &
    buffer_send_vector_gpu,buffer_recv_vector_gpu, &
    b_buffer_send_vector_gpu,b_buffer_recv_vector_gpu, &
    request_send_recv_vector_gpu,b_request_send_recv_vector_gpu

  implicit none

  logical,intent(in) :: compute_b_wavefield

  ! local parameters
  integer:: iphase

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

  ! check
  if (PML_BOUNDARY_CONDITIONS ) &
    call exit_MPI(myrank,'PML conditions not yet implemented for routine compute_forces_viscoelastic_GPU()')

  ! distinguishes two runs: for points on MPI interfaces, and points within the partitions
  do iphase = 1,2

    ! elastic term
    ! contains both forward SIM_TYPE==1 and backward SIM_TYPE==3 simulations
    call compute_forces_viscoelastic_cuda(Mesh_pointer, iphase, deltat, &
                                          nspec_outer_elastic, &
                                          nspec_inner_elastic, &
                                          any_anisotropy,ATTENUATION_VISCOELASTIC, &
                                          compute_wavefield_1,compute_wavefield_2)


    ! while inner elements compute "Kernel_2", we wait for MPI to
    ! finish and transfer the boundary terms to the device asynchronously
    if (NPROC > 1) then
      if (iphase == 2) then
        !daniel: todo - this avoids calling the Fortran vector send from CUDA routine
        if (compute_wavefield_1) then
          ! wait for asynchronous copy to finish
          call sync_copy_from_device(Mesh_pointer,iphase,buffer_send_vector_gpu)

          ! sends MPI buffers
          call assemble_MPI_vector_send_cuda(NPROC, &
                                             buffer_send_vector_gpu,buffer_recv_vector_gpu, &
                                             ninterface,max_nibool_interfaces_ext_mesh, &
                                             nibool_interfaces_ext_mesh, &
                                             my_neighbors, &
                                             request_send_recv_vector_gpu,ninterface_elastic,inum_interfaces_elastic)

          ! transfers MPI buffers onto GPU
          call transfer_boundary_to_device(NPROC,Mesh_pointer,buffer_recv_vector_gpu, &
                                           ninterface,max_nibool_interfaces_ext_mesh, &
                                           request_send_recv_vector_gpu,ninterface_elastic,inum_interfaces_elastic)
        endif
      endif ! inner elements
    endif

    ! computes additional contributions
    if (iphase == 1) then
      ! adds elastic absorbing boundary term to acceleration (Stacey conditions)
      if (STACEY_ABSORBING_CONDITIONS) then
        call compute_stacey_viscoelastic_GPU(iphase,compute_wavefield_1,compute_wavefield_2)
      endif

      ! acoustic coupling
      if (coupled_acoustic_elastic) then
        if (compute_wavefield_1) &
          call compute_coupling_el_ac_cuda(Mesh_pointer,iphase,num_fluid_solid_edges,1) ! 1 == forward/adjoint
        if (compute_wavefield_2) &
          call compute_coupling_el_ac_cuda(Mesh_pointer,iphase,num_fluid_solid_edges,3) ! 3 == backward
      endif

      ! poroelastic coupling
      if (coupled_elastic_poro) then
        call stop_the_code('Error GPU simulation: poroelastic coupling not implemented yet')
      endif

      ! adds source term (single-force/moment-tensor solution)
      call compute_add_sources_viscoelastic_GPU(iphase,compute_wavefield_1,compute_wavefield_2)
    endif

    ! assemble all the contributions between slices using MPI
    if (NPROC > 1) then
      if (iphase == 1) then
        ! sends accel values to corresponding MPI interface neighbors
        if (compute_wavefield_1) then
          ! transfers boundary region to host asynchronously. The
          ! MPI-send is done from within compute_forces_viscoelastic_cuda,
          ! once the inner element kernels are launched, and the
          ! memcpy has finished. see compute_forces_viscoelastic_cuda: ~ line 1655
          call transfer_boundary_from_device_a(Mesh_pointer)
        endif
        ! adjoint simulations
        if (compute_wavefield_2) then
           call transfer_boun_accel_from_device(Mesh_pointer, &
                        b_buffer_send_vector_gpu, &
                        3) ! -- 3 == adjoint b_accel
           call assemble_MPI_vector_send_cuda(NPROC, &
                        b_buffer_send_vector_gpu,b_buffer_recv_vector_gpu, &
                        ninterface,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh, &
                        my_neighbors, &
                        b_request_send_recv_vector_gpu,ninterface_elastic,inum_interfaces_elastic)
        endif !adjoint

      else
        ! waits for send/receive requests to be completed and assembles values
        if (compute_wavefield_1) then
          call assemble_MPI_vector_write_cuda(NPROC,Mesh_pointer, &
                                              buffer_recv_vector_gpu,ninterface, &
                                              max_nibool_interfaces_ext_mesh, &
                                              nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                              request_send_recv_vector_gpu, &
                                              1,ninterface_elastic,inum_interfaces_elastic)
        endif
        ! adjoint simulations
        if (compute_wavefield_2) then
          call assemble_MPI_vector_write_cuda(NPROC,Mesh_pointer, &
                                              b_buffer_recv_vector_gpu,ninterface, &
                                              max_nibool_interfaces_ext_mesh, &
                                              nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                              b_request_send_recv_vector_gpu, &
                                              3,ninterface_elastic,inum_interfaces_elastic)
        endif !adjoint
      endif
    endif ! NPROC > 1

  enddo

  ! multiplies with inverse of mass matrix (note: rmass has been inverted already)
  call kernel_3_a_cuda(Mesh_pointer,deltatover2,b_deltatover2,compute_wavefield_1,compute_wavefield_2)

  end subroutine compute_forces_viscoelastic_GPU

!
!---------------------------------------------------------------------------------------------
!

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic_GPU(iphase,compute_wavefield_1,compute_wavefield_2)

  use constants, only: CUSTOM_REAL,NGLLX,NDIM

  use specfem_par, only: nspec_bottom,nspec_left,nspec_top,nspec_right, &
                         b_absorb_elastic_left,b_absorb_elastic_right, &
                         b_absorb_elastic_bottom, b_absorb_elastic_top, &
                         SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                         num_abs_boundary_faces,NO_BACKWARD_RECONSTRUCTION,UNDO_ATTENUATION_AND_OR_PML

  use specfem_par_gpu, only: Mesh_pointer

  implicit none

  ! communication overlap
  integer,intent(in) :: iphase
  logical,intent(in) :: compute_wavefield_1,compute_wavefield_2

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,nspec_bottom) :: b_absorb_elastic_bottom_slice
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,nspec_left) :: b_absorb_elastic_left_slice
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,nspec_right) :: b_absorb_elastic_right_slice
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,nspec_top) :: b_absorb_elastic_top_slice

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return
  if (NO_BACKWARD_RECONSTRUCTION) return

  if (SIMULATION_TYPE == 3) then
    ! gets absorbing contribution buffers
    b_absorb_elastic_bottom_slice(:,:,:) = b_absorb_elastic_bottom(:,:,:,NSTEP-it+1)
    b_absorb_elastic_left_slice(:,:,:) = b_absorb_elastic_left(:,:,:,NSTEP-it+1)
    b_absorb_elastic_right_slice(:,:,:) = b_absorb_elastic_right(:,:,:,NSTEP-it+1)
    b_absorb_elastic_top_slice(:,:,:) = b_absorb_elastic_top(:,:,:,NSTEP-it+1)
  endif

  call compute_stacey_viscoelastic_cuda(Mesh_pointer,iphase, &
                                        b_absorb_elastic_left_slice,b_absorb_elastic_right_slice, &
                                        b_absorb_elastic_top_slice,b_absorb_elastic_bottom_slice, &
                                        compute_wavefield_1,compute_wavefield_2, &
                                        UNDO_ATTENUATION_AND_OR_PML)

  ! adjoint simulations: stores absorbed wavefield part
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
    ! writes out absorbing boundary value only when second phase is running
    b_absorb_elastic_bottom(:,:,:,it) = b_absorb_elastic_bottom_slice(:,:,:)
    b_absorb_elastic_right(:,:,:,it) = b_absorb_elastic_right_slice(:,:,:)
    b_absorb_elastic_top(:,:,:,it) = b_absorb_elastic_top_slice(:,:,:)
    b_absorb_elastic_left(:,:,:,it) = b_absorb_elastic_left_slice(:,:,:)
  endif

  end subroutine compute_stacey_viscoelastic_GPU

!
!---------------------------------------------------------------------------------------------
!

! for elastic solver on GPU

  subroutine compute_add_sources_viscoelastic_GPU(iphase,compute_wavefield_1,compute_wavefield_2)

  use specfem_par, only: NSTEP,SIMULATION_TYPE,nadj_rec_local,it

  use specfem_par_gpu, only: Mesh_pointer

  implicit none

  integer,intent(in) :: iphase
  logical,intent(in) :: compute_wavefield_1,compute_wavefield_2

  ! local parameters
  integer :: it_tmp

  ! forward simulations
  if (SIMULATION_TYPE == 1) call compute_add_sources_el_cuda(Mesh_pointer,iphase,it)

  ! adjoint simulations
  ! time step index
  it_tmp = NSTEP - it + 1

  ! adds adjoint sources
  if (SIMULATION_TYPE /= 1 .and. nadj_rec_local > 0 .and. compute_wavefield_1) &
    call add_sources_el_sim_type_2_or_3(Mesh_pointer, iphase, it_tmp, nadj_rec_local, NSTEP)

  ! kernel simulations w/ backward/reconstructed sources
  if (compute_wavefield_2) call compute_add_sources_el_s3_cuda(Mesh_pointer, iphase, it_tmp)

  end subroutine compute_add_sources_viscoelastic_GPU
