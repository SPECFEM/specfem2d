subroutine compute_forces_acoustic_GPU(STACEY_BOUNDARY_CONDITIONS,nspec,nelemabs,ninterface, &
                           ninterface_acoustic,inum_interfaces_acoustic)

  use specfem_par


  implicit none

   !Parametres d'entree

  integer :: nspec, nelemabs, ninterface, ninterface_acoustic
  logical :: STACEY_BOUNDARY_CONDITIONS
integer, dimension(ninterface), intent(in)  :: inum_interfaces_acoustic


  


  ! local parameters
  integer:: iphase,i,j
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
      call compute_stacey_acoustic_GPU(phase_is_inner,nelemabs,&
                            SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                            Mesh_pointer)
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
          stop 'not implemented yet'
    endif


    ! sources
    call compute_add_sources_acoustic_GPU(phase_is_inner, &
                                  NSOURCES,it,&
                                  SIMULATION_TYPE,NSTEP, &
                                  nadj_rec_local, &
                                  Mesh_pointer)

! assemble all the contributions between slices using MPI

    if(ninterface_acoustic > 0) then

    if( phase_is_inner .eqv. .false. ) then

      ! sends potential_dot_dot_acoustic values to corresponding MPI interface neighbors (non-blocking)

     call transfer_boun_pot_from_device(Mesh_pointer, &
                                         buffer_send_scalar_ext_mesh, &
                                         1) ! <-- 1 == fwd accel



!buffer_recv_scalar_ext_mesh(:,:)=0

      call assemble_MPI_scalar_send_cuda(NPROC, &
                        buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                        ninterface,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,&
                        my_neighbours, &
                        tab_requests_send_recv_scalar,ninterface_acoustic,inum_interfaces_acoustic)




      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) then
        call transfer_boun_pot_from_device(Mesh_pointer, &
                                           buffer_send_scalar_ext_mesh,&
                                           3) ! <-- 3 == adjoint b_accel

        call assemble_MPI_scalar_send_cuda(NPROC, &
                          buffer_send_scalar_ext_mesh,b_buffer_recv_scalar_ext_mesh, &
                          ninterface,max_nibool_interfaces_ext_mesh, &
                          nibool_interfaces_ext_mesh,&
                          my_neighbours, &
                          b_tab_requests_send_recv_scalar,ninterface_acoustic,inum_interfaces_acoustic)

      endif

    else

      ! waits for send/receive requests to be completed and assembles values
      call assemble_MPI_scalar_write_cuda(NPROC,NGLOB_AB, &
                        Mesh_pointer,&
                        buffer_recv_scalar_ext_mesh, &
                        ninterface, &
                        max_nibool_interfaces_ext_mesh, &
                        tab_requests_send_recv_scalar, &
                        1,ninterface_acoustic,inum_interfaces_acoustic)



      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) then
        call assemble_MPI_scalar_write_cuda(NPROC,NGLOB_AB, &
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
  call kernel_3_a_acoustic_cuda(Mesh_pointer)
 
    
! corrector:
! updates the chi_dot term which requires chi_dot_dot(t+delta)
  call kernel_3_b_acoustic_cuda(Mesh_pointer,deltatover2f,b_deltatover2f)

! enforces free surface (zeroes potentials at free surface)

 call acoustic_enforce_free_surf_cuda(Mesh_pointer)


end subroutine compute_forces_acoustic_GPU
!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! for acoustic solver on GPU

  subroutine compute_stacey_acoustic_GPU(phase_is_inner,num_abs_boundary_faces,&
                            SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                            Mesh_pointer)

  use constants
  use specfem_par, only : nspec_bottom,nspec_left,nspec_top,nspec_right,b_absorb_acoustic_left,b_absorb_acoustic_right,&
                           b_absorb_acoustic_bottom, b_absorb_acoustic_top,myrank

  implicit none

! potentials

! communication overlap
  logical :: phase_is_inner

! absorbing boundary surface
  integer :: num_abs_boundary_faces

! adjoint simulations
  integer:: SIMULATION_TYPE
  integer:: NSTEP,it
  logical:: SAVE_FORWARD
  real(kind=CUSTOM_REAL),dimension(NGLLX,nspec_bottom) :: b_absorb_potential_bottom_slice
  real(kind=CUSTOM_REAL),dimension(NGLLX,nspec_left) :: b_absorb_potential_left_slice
  real(kind=CUSTOM_REAL),dimension(NGLLX,nspec_right) :: b_absorb_potential_right_slice
  real(kind=CUSTOM_REAL),dimension(NGLLX,nspec_top) :: b_absorb_potential_top_slice


  ! GPU_MODE variables
  integer(kind=8) :: Mesh_pointer

  ! checks if anything to do
  if( num_abs_boundary_faces == 0 ) return
  

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

  subroutine compute_add_sources_acoustic_GPU(phase_is_inner, &
                                  NSOURCES,it,&
                                  SIMULATION_TYPE,NSTEP, &
                                  nadj_rec_local, &
                                  Mesh_pointer)

  use constants
  use specfem_par,only: nsources_local
                        
  implicit none


! communication overlap
  logical :: phase_is_inner

! source
  integer :: NSOURCES,it



!adjoint simulations
  integer:: SIMULATION_TYPE,NSTEP
  integer(kind=8) :: Mesh_pointer
  integer:: nadj_rec_local


! forward simulations
  if (SIMULATION_TYPE == 1 .and. nsources_local > 0) then
      call compute_add_sources_ac_cuda(Mesh_pointer,phase_is_inner,NSOURCES,it)
  endif

! adjoint simulations
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! adds adjoint source in this partitions
    if( nadj_rec_local > 0 ) then
      if( it < NSTEP ) then
        ! receivers act as sources
        ! on GPU
        call add_sources_ac_sim_2_or_3_cuda(Mesh_pointer,phase_is_inner, NSTEP -it + 1, nadj_rec_local,NSTEP)
      endif ! it
    endif ! nadj_rec_local > 0
  endif

! adjoint simulations
  if (SIMULATION_TYPE == 3 .and. nsources_local > 0) then
      call compute_add_sources_ac_s3_cuda(Mesh_pointer,phase_is_inner,NSOURCES,NSTEP -it + 1)
  endif

  end subroutine compute_add_sources_acoustic_GPU
