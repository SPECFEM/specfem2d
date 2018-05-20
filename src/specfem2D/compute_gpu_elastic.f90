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

  subroutine compute_forces_viscoelastic_GPU()

  use constants, only: CUSTOM_REAL,NGLLX,NDIM

  use specfem_par, only: myrank,NPROC,ninterface,max_nibool_interfaces_ext_mesh,nibool_interfaces_ext_mesh, &
    my_neighbors,ninterface_elastic,inum_interfaces_elastic,ibool_interfaces_ext_mesh, &
    num_fluid_solid_edges,nspec_bottom,nspec_left,nspec_right,nspec_top, &
    STACEY_ABSORBING_CONDITIONS,PML_BOUNDARY_CONDITIONS,any_poroelastic,any_acoustic,SIMULATION_TYPE,ATTENUATION_VISCOELASTIC

  use specfem_par, only: nspec_outer_elastic,nspec_inner_elastic

  use specfem_par_gpu, only: Mesh_pointer,ANY_ANISOTROPY,deltatf,deltatover2f,b_deltatover2f, &
    buffer_send_vector_gpu,buffer_recv_vector_gpu, &
    b_buffer_send_vector_gpu,b_buffer_recv_vector_gpu, &
    request_send_recv_vector_gpu,b_request_send_recv_vector_gpu

  implicit none

  ! local parameters
  integer:: iphase
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,nspec_bottom) :: b_absorb_elastic_bottom_slice
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,nspec_left) :: b_absorb_elastic_left_slice
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,nspec_right) :: b_absorb_elastic_right_slice
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,nspec_top) :: b_absorb_elastic_top_slice

  ! check
  if (PML_BOUNDARY_CONDITIONS ) &
    call exit_MPI(myrank,'PML conditions not yet implemented for routine compute_forces_viscoelastic_GPU()')

  ! distinguishes two runs: for points on MPI interfaces, and points within the partitions
  do iphase = 1,2

    ! elastic term
    ! contains both forward SIM_TYPE==1 and backward SIM_TYPE==3 simulations
    call compute_forces_viscoelastic_cuda(Mesh_pointer, iphase, deltatf, &
                                          nspec_outer_elastic, &
                                          nspec_inner_elastic, &
                                          ANY_ANISOTROPY,ATTENUATION_VISCOELASTIC)


    ! while inner elements compute "Kernel_2", we wait for MPI to
    ! finish and transfer the boundary terms to the device asynchronously
    if (NPROC > 1) then
      if (iphase == 2) then
        !daniel: todo - this avoids calling the Fortran vector send from CUDA routine
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
      endif ! inner elements
    endif

    ! computes additional contributions
    if (iphase == 1) then
      ! adds elastic absorbing boundary term to acceleration (Stacey conditions)
      if (STACEY_ABSORBING_CONDITIONS) then
        call compute_stacey_viscoelastic_GPU(iphase,b_absorb_elastic_bottom_slice,b_absorb_elastic_left_slice, &
                                             b_absorb_elastic_right_slice,b_absorb_elastic_top_slice)
      endif

      ! acoustic coupling
      if (any_acoustic) then
        if (num_fluid_solid_edges > 0) then
          call compute_coupling_el_ac_cuda(Mesh_pointer,iphase, &
                                           num_fluid_solid_edges)
        endif
      endif

      ! poroelastic coupling
      ! poroelastic coupling
      if (any_poroelastic) then
        call stop_the_code('Error GPU simulation: poroelastic coupling not implemented yet')
      endif

      ! adds source term (single-force/moment-tensor solution)
      call compute_add_sources_viscoelastic_GPU(iphase)
    endif

    ! assemble all the contributions between slices using MPI
    if (NPROC > 1) then
      if (iphase == 1) then
        ! sends accel values to corresponding MPI interface neighbors

        ! transfers boundary region to host asynchronously. The
        ! MPI-send is done from within compute_forces_viscoelastic_cuda,
        ! once the inner element kernels are launched, and the
        ! memcpy has finished. see compute_forces_viscoelastic_cuda: ~ line 1655
        call transfer_boundary_from_device_a(Mesh_pointer)

        ! adjoint simulations
        if (SIMULATION_TYPE == 3) then
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
        call assemble_MPI_vector_write_cuda(NPROC,Mesh_pointer, &
                                            buffer_recv_vector_gpu,ninterface, &
                                            max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            request_send_recv_vector_gpu, &
                                            1,ninterface_elastic,inum_interfaces_elastic)

        ! adjoint simulations
        if (SIMULATION_TYPE == 3) then
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
  call kernel_3_a_cuda(Mesh_pointer,deltatover2f,b_deltatover2f)

  end subroutine compute_forces_viscoelastic_GPU

!
!---------------------------------------------------------------------------------------------
!

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic_GPU(iphase,b_absorb_elastic_bottom_slice,b_absorb_elastic_left_slice, &
                                             b_absorb_elastic_right_slice, b_absorb_elastic_top_slice)

  use constants, only: CUSTOM_REAL,NGLLX,NDIM

  use specfem_par, only: nspec_bottom,nspec_left,nspec_top,nspec_right,b_absorb_elastic_left,b_absorb_elastic_right, &
                          b_absorb_elastic_bottom, b_absorb_elastic_top,SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                          nelemabs

  use specfem_par_gpu, only: Mesh_pointer

  implicit none

  ! communication overlap
  integer :: iphase

  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,nspec_bottom) :: b_absorb_elastic_bottom_slice
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,nspec_left) :: b_absorb_elastic_left_slice
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,nspec_right) :: b_absorb_elastic_right_slice
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,nspec_top) :: b_absorb_elastic_top_slice

  ! checks if anything to do
  if (nelemabs == 0) return

  if (SIMULATION_TYPE == 3) then
    ! gets absorbing contribution buffers
    b_absorb_elastic_bottom_slice(1,:,:) = b_absorb_elastic_bottom(1,:,:,NSTEP-it+1)
    b_absorb_elastic_left_slice(1,:,:) = b_absorb_elastic_left(1,:,:,NSTEP-it+1)
    b_absorb_elastic_right_slice(1,:,:) = b_absorb_elastic_right(1,:,:,NSTEP-it+1)
    b_absorb_elastic_top_slice(1,:,:) = b_absorb_elastic_top(1,:,:,NSTEP-it+1)

    b_absorb_elastic_bottom_slice(2,:,:) = b_absorb_elastic_bottom(2,:,:,NSTEP-it+1)
    b_absorb_elastic_left_slice(2,:,:) = b_absorb_elastic_left(2,:,:,NSTEP-it+1)
    b_absorb_elastic_right_slice(2,:,:) = b_absorb_elastic_right(2,:,:,NSTEP-it+1)
    b_absorb_elastic_top_slice(2,:,:) = b_absorb_elastic_top(2,:,:,NSTEP-it+1)
  endif

  call compute_stacey_viscoelastic_cuda(Mesh_pointer,iphase,b_absorb_elastic_left_slice, &
                   b_absorb_elastic_right_slice,b_absorb_elastic_top_slice,b_absorb_elastic_bottom_slice)

  ! adjoint simulations: stores absorbed wavefield part
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
    ! writes out absorbing boundary value only when second phase is running
    b_absorb_elastic_bottom(1,:,:,it) = b_absorb_elastic_bottom_slice(1,:,:)
    b_absorb_elastic_right(1,:,:,it) = b_absorb_elastic_right_slice(1,:,:)
    b_absorb_elastic_top(1,:,:,it) = b_absorb_elastic_top_slice(1,:,:)
    b_absorb_elastic_left(1,:,:,it) = b_absorb_elastic_left_slice(1,:,:)

    b_absorb_elastic_bottom(2,:,:,it) = b_absorb_elastic_bottom_slice(2,:,:)
    b_absorb_elastic_right(2,:,:,it) = b_absorb_elastic_right_slice(2,:,:)
    b_absorb_elastic_top(2,:,:,it) = b_absorb_elastic_top_slice(2,:,:)
    b_absorb_elastic_left(2,:,:,it) = b_absorb_elastic_left_slice(2,:,:)
  endif

  end subroutine compute_stacey_viscoelastic_GPU

!
!---------------------------------------------------------------------------------------------
!

! for elastic solver on GPU

  subroutine compute_add_sources_viscoelastic_GPU(iphase)

  use specfem_par, only: NSTEP,SIMULATION_TYPE,nadj_rec_local,it

  use specfem_par_gpu, only: Mesh_pointer

  implicit none

  integer,intent(in) :: iphase

  ! forward simulations
  if (SIMULATION_TYPE == 1) call compute_add_sources_el_cuda(Mesh_pointer,iphase,it)

  ! adjoint simulations
  if (SIMULATION_TYPE == 3 .and. nadj_rec_local > 0 .and. it < NSTEP) then
    call add_sources_el_sim_type_2_or_3(Mesh_pointer,iphase, NSTEP -it + 1, nadj_rec_local,NSTEP)
  endif

  ! adjoint simulations
  if (SIMULATION_TYPE == 3) call compute_add_sources_el_s3_cuda(Mesh_pointer,iphase,NSTEP -it + 1)

  end subroutine compute_add_sources_viscoelastic_GPU
