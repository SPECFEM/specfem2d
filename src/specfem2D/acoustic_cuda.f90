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



subroutine compute_forces_acoustic_GPU()

  use specfem_par, only : Mesh_pointer,nspec_outer_acoustic, nspec_inner_acoustic,num_fluid_solid_edges, &
                          buffer_send_scalar_ext_mesh,NPROC,buffer_recv_scalar_ext_mesh,deltatover2f,b_deltatover2f, &
                          ninterface,max_nibool_interfaces_ext_mesh,nibool_interfaces_ext_mesh, &
                          my_neighbours,tab_requests_send_recv_scalar,ninterface_acoustic,inum_interfaces_acoustic, &
                          b_buffer_send_scalar_ext_mesh,b_buffer_recv_scalar_ext_mesh,b_tab_requests_send_recv_scalar, &
                          nelem_acoustic_surface,STACEY_BOUNDARY_CONDITIONS,any_elastic,any_poroelastic,SIMULATION_TYPE

  implicit none


  ! local parameters
  integer:: iphase
  logical:: phase_is_inner


   if(nelem_acoustic_surface > 0) then

   !  enforces free surface (zeroes potentials at free surface)
    call acoustic_enforce_free_surf_cuda(Mesh_pointer)

  endif


  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase=1,2


    !first for points on MPI interfaces, thus outer elements
    if( iphase == 1 ) then
      phase_is_inner = .false.
    else
      phase_is_inner = .true.
    endif



    ! acoustic pressure term
    ! includes code for SIMULATION_TYPE==3
    call compute_forces_acoustic_cuda(Mesh_pointer, iphase, &
                                      nspec_outer_acoustic, nspec_inner_acoustic)




    ! ! Stacey absorbing boundary conditions
    if(STACEY_BOUNDARY_CONDITIONS) then
      call compute_stacey_acoustic_GPU(phase_is_inner)
    endif



    ! elastic coupling
    if(any_elastic ) then
      if( num_fluid_solid_edges > 0 ) then
        ! on GPU
        call compute_coupling_ac_el_cuda(Mesh_pointer,phase_is_inner, &
                                         num_fluid_solid_edges)
      endif
    endif

! poroelastic coupling
    if(any_poroelastic )  then
          stop 'poroelastic not implemented yet in GPU mode'
    endif


    ! sources
    call compute_add_sources_acoustic_GPU(phase_is_inner)



! assemble all the contributions between slices using MPI

    if(ninterface_acoustic > 0) then

    if( phase_is_inner .eqv. .false. ) then

      ! sends potential_dot_dot_acoustic values to corresponding MPI interface neighbors (non-blocking)

     call transfer_boun_pot_from_device(Mesh_pointer, &
                                         buffer_send_scalar_ext_mesh, &
                                         1) ! <-- 1 == fwd accel

      call assemble_MPI_scalar_send_cuda(NPROC, &
                        buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                        ninterface,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,&
                        my_neighbours, &
                        tab_requests_send_recv_scalar,ninterface_acoustic,inum_interfaces_acoustic)




      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) then
        call transfer_boun_pot_from_device(Mesh_pointer, &
                                           b_buffer_send_scalar_ext_mesh,&
                                           3) ! <-- 3 == adjoint b_accel

        call assemble_MPI_scalar_send_cuda(NPROC, &
                          b_buffer_send_scalar_ext_mesh,b_buffer_recv_scalar_ext_mesh, &
                          ninterface,max_nibool_interfaces_ext_mesh, &
                          nibool_interfaces_ext_mesh,&
                          my_neighbours, &
                          b_tab_requests_send_recv_scalar,ninterface_acoustic,inum_interfaces_acoustic)

      endif

    else


      ! waits for send/receive requests to be completed and assembles values
      call assemble_MPI_scalar_write_cuda(NPROC, &
                        Mesh_pointer,&
                        buffer_recv_scalar_ext_mesh, &
                        ninterface, &
                        max_nibool_interfaces_ext_mesh, &
                        tab_requests_send_recv_scalar, &
                        1,ninterface_acoustic,inum_interfaces_acoustic)



      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) then
        call assemble_MPI_scalar_write_cuda(NPROC, &
                        Mesh_pointer, &
                        b_buffer_recv_scalar_ext_mesh, &
                        ninterface, &
                        max_nibool_interfaces_ext_mesh, &
                        b_tab_requests_send_recv_scalar, &
                        3,ninterface_acoustic,inum_interfaces_acoustic)
      endif
    endif !phase_is_inner

  endif !interface_acoustic


  enddo


! divides pressure with mass matrix
! corrector:
! updates the chi_dot term which requires chi_dot_dot(t+delta)
  call kernel_3_acoustic_cuda(Mesh_pointer,deltatover2f,b_deltatover2f)

! enforces free surface (zeroes potentials at free surface)
 call acoustic_enforce_free_surf_cuda(Mesh_pointer)


end subroutine compute_forces_acoustic_GPU
!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! for acoustic solver on GPU

  subroutine compute_stacey_acoustic_GPU(phase_is_inner)

  use constants
  use specfem_par, only : nspec_bottom,nspec_left,nspec_top,nspec_right,b_absorb_acoustic_left,b_absorb_acoustic_right,&
                          b_absorb_acoustic_bottom, b_absorb_acoustic_top,it,NSTEP,SIMULATION_TYPE,SAVE_FORWARD,&
                          nelemabs,Mesh_pointer

  implicit none

! communication overlap
  logical :: phase_is_inner

! adjoint simulations
  real(kind=CUSTOM_REAL),dimension(NGLLX,nspec_bottom) :: b_absorb_potential_bottom_slice
  real(kind=CUSTOM_REAL),dimension(NGLLX,nspec_left) :: b_absorb_potential_left_slice
  real(kind=CUSTOM_REAL),dimension(NGLLX,nspec_right) :: b_absorb_potential_right_slice
  real(kind=CUSTOM_REAL),dimension(NGLLX,nspec_top) :: b_absorb_potential_top_slice

  ! checks if anything to do
  if( nelemabs == 0 ) return


if( SIMULATION_TYPE == 3 ) then
    if( phase_is_inner .eqv. .false. ) then
    b_absorb_potential_left_slice(:,:)=b_absorb_acoustic_left(:,:,NSTEP-it+1)
    b_absorb_potential_right_slice(:,:)=b_absorb_acoustic_right(:,:,NSTEP-it+1)
    b_absorb_potential_top_slice(:,:)=b_absorb_acoustic_top(:,:,NSTEP-it+1)
    b_absorb_potential_bottom_slice(:,:)=b_absorb_acoustic_bottom(:,:,NSTEP-it+1)
    endif
endif


  ! absorbs absorbing-boundary surface using Sommerfeld condition (vanishing field in the outer-space)
  call compute_stacey_acoustic_cuda(Mesh_pointer, phase_is_inner,b_absorb_potential_left_slice,b_absorb_potential_right_slice,&
                                     b_absorb_potential_top_slice,b_absorb_potential_bottom_slice)

  ! adjoint simulations: stores absorbed wavefield part
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD ) then
    ! writes out absorbing boundary value only when second phase is running
    if( phase_is_inner .eqv. .true. ) then

    b_absorb_acoustic_bottom(:,:,it) = b_absorb_potential_bottom_slice(:,:)
    b_absorb_acoustic_right(:,:,it) = b_absorb_potential_right_slice(:,:)
    b_absorb_acoustic_top(:,:,it) = b_absorb_potential_top_slice(:,:)
    b_absorb_acoustic_left(:,:,it) = b_absorb_potential_left_slice(:,:)

    endif
  endif

  end subroutine compute_stacey_acoustic_GPU


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine compute_add_sources_acoustic_GPU(phase_is_inner)

  use constants
  use specfem_par,only: it,SIMULATION_TYPE,NSTEP, nadj_rec_local,Mesh_pointer

  implicit none

! communication overlap
  logical :: phase_is_inner

! forward simulations
  if (SIMULATION_TYPE == 1 ) call compute_add_sources_ac_cuda(Mesh_pointer,phase_is_inner,it)

! adjoint simulations
  if ( SIMULATION_TYPE == 3 .and. nadj_rec_local > 0 .and. it < NSTEP ) then
        call add_sources_ac_sim_2_or_3_cuda(Mesh_pointer,phase_is_inner, NSTEP -it + 1, nadj_rec_local,NSTEP)
  endif

! adjoint simulations
  if (SIMULATION_TYPE == 3) call compute_add_sources_ac_s3_cuda(Mesh_pointer,phase_is_inner,NSTEP -it + 1)

  end subroutine compute_add_sources_acoustic_GPU
