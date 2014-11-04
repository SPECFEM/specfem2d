
! elastic solver

subroutine compute_forces_elastic_GPU(STACEY_BOUNDARY_CONDITIONS,ninterface,nelemabs,ninterface_elastic,inum_interfaces_elastic)

  use specfem_par

  implicit none


  integer :: ninterface, nelemabs,i,j, ninterface_elastic
  integer, dimension(ninterface), intent(in)  :: inum_interfaces_elastic
  integer:: iphase
  logical:: phase_is_inner,STACEY_BOUNDARY_CONDITIONS
  real(kind=CUSTOM_REAL),dimension(2,NGLLX,nspec_bottom) :: b_absorb_elastic_bottom_slice
  real(kind=CUSTOM_REAL),dimension(2,NGLLX,nspec_left) :: b_absorb_elastic_left_slice
  real(kind=CUSTOM_REAL),dimension(2,NGLLX,nspec_right) :: b_absorb_elastic_right_slice
  real(kind=CUSTOM_REAL),dimension(2,NGLLX,nspec_top) :: b_absorb_elastic_top_slice



  ! check
  if( PML_BOUNDARY_CONDITIONS ) &
    call exit_MPI(myrank,'PML conditions not yet implemented for routine compute_forces_viscoelastic_GPU()')

  ! distinguishes two runs: for points on MPI interfaces, and points within the partitions
  do iphase=1,2

    !first for points on MPI interfaces
    if( iphase == 1 ) then
      phase_is_inner = .false.
    else
      phase_is_inner = .true.
    endif

    ! elastic term
    ! contains both forward SIM_TYPE==1 and backward SIM_TYPE==3 simulations
    call compute_forces_viscoelastic_cuda(Mesh_pointer, iphase, deltatf, &
                                          nspec_outer_elastic, &
                                          nspec_inner_elastic, &
                                          ANY_ANISOTROPY)


    ! while inner elements compute "Kernel_2", we wait for MPI to
    ! finish and transfer the boundary terms to the device asynchronously
    if(phase_is_inner .eqv. .true.) then
      !daniel: todo - this avoids calling the fortran vector send from CUDA routine
      ! wait for asynchronous copy to finish

    if(.not. CUDA_AWARE_MPI) then  

     call sync_copy_from_device(Mesh_pointer,iphase,buffer_send_vector_ext_mesh)


      ! sends mpi buffers
      call assemble_MPI_vector_send_cuda(NPROC, &
                  buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                  ninterface,max_nibool_interfaces_ext_mesh, &
                  nibool_interfaces_ext_mesh,&
                  my_neighbours, &
                  tab_requests_send_recv_vector,ninterface_elastic,inum_interfaces_elastic)




      ! transfers mpi buffers onto GPU
      call transfer_boundary_to_device(NPROC,Mesh_pointer,buffer_recv_vector_ext_mesh, &
                  ninterface,max_nibool_interfaces_ext_mesh, &
                  tab_requests_send_recv_vector,ninterface_elastic,inum_interfaces_elastic)

     endif ! not CUDA_AWARE_MPI

    endif ! inner elements

    ! adds elastic absorbing boundary term to acceleration (Stacey conditions)

    if( STACEY_BOUNDARY_CONDITIONS ) then
      call compute_stacey_viscoelastic_GPU(phase_is_inner,nelemabs, &
                                           SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                                           Mesh_pointer,b_absorb_elastic_bottom_slice,b_absorb_elastic_left_slice, &
                                           b_absorb_elastic_right_slice,b_absorb_elastic_top_slice)
    endif

    ! acoustic coupling
    if( any_acoustic ) then
      if( num_fluid_solid_edges > 0 ) then
        call compute_coupling_el_ac_cuda(Mesh_pointer,phase_is_inner, &
                                         num_fluid_solid_edges)
      endif
    endif

    ! poroelastic coupling
    ! poroelastic coupling
    if(any_poroelastic )  then
          stop 'not implemented yet'
    endif


    ! adds source term (single-force/moment-tensor solution)


    call compute_add_sources_viscoelastic_GPU(phase_is_inner, &
                                              NSOURCES,it,&
                                              SIMULATION_TYPE,NSTEP, &
                                              nadj_rec_local, &
                                              Mesh_pointer)


    ! assemble all the contributions between slices using MPI
    if( phase_is_inner .eqv. .false. ) then
  
      if(.not. CUDA_AWARE_MPI) then 
        ! sends accel values to corresponding MPI interface neighbors

        ! transfers boundary region to host asynchronously. The
        ! MPI-send is done from within compute_forces_viscoelastic_cuda,
        ! once the inner element kernels are launched, and the
        ! memcpy has finished. see compute_forces_viscoelastic_cuda: ~ line 1655



        call transfer_boundary_from_device_a(Mesh_pointer)

        ! adjoint simulations
        if( SIMULATION_TYPE == 3 ) then
           call transfer_boun_accel_from_device(Mesh_pointer,&
                        b_buffer_send_vector_ext_mesh,&
                        3) ! <-- 3 == adjoint b_accel
           call assemble_MPI_vector_send_cuda(NPROC, &
                        b_buffer_send_vector_ext_mesh,b_buffer_recv_vector_ext_mesh, &
                        ninterface,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,&
                        my_neighbours, &
                        b_tab_requests_send_recv_vector,ninterface_elastic,inum_interfaces_elastic)
        endif !adjoint

      


       else

      
         call prepare_boundary_on_device(Mesh_pointer)

      
         call transfer_accel_aware(Mesh_pointer,NPROC,max_nibool_interfaces_ext_mesh, &
                  nibool_interfaces_ext_mesh,&
                  my_neighbours, &
                  tab_requests_send_recv_vector_c,inum_interfaces_elastic)
 
      

       endif !CUDA_AWARE_MPI

    else

       if(.not. CUDA_AWARE_MPI) then 
!       waits for send/receive requests to be completed and assembles values
         call assemble_MPI_vector_write_cuda(NPROC,NGLOB_AB, Mesh_pointer,&
                      buffer_recv_vector_ext_mesh,ninterface,&
                      max_nibool_interfaces_ext_mesh, &
                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                      tab_requests_send_recv_vector, &
                      1,ninterface_elastic,inum_interfaces_elastic)
         ! adjoint simulations
         if( SIMULATION_TYPE == 3 ) then
         call assemble_MPI_vector_write_cuda(NPROC,NGLOB_AB, Mesh_pointer,&
                              b_buffer_recv_vector_ext_mesh,ninterface,&
                              max_nibool_interfaces_ext_mesh, &
                              nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                              b_tab_requests_send_recv_vector, &
                              3,ninterface_elastic,inum_interfaces_elastic)
         endif !adjoint

       else

        call receive_accel_aware(tab_requests_send_recv_vector_c,ninterface_elastic,NPROC,inum_interfaces_elastic)

        call transfer_asmbl_accel_to_device(Mesh_pointer, &
                                        buffer_recv_vector_ext_mesh, &
                                        max_nibool_interfaces_ext_mesh, &
                                        nibool_interfaces_ext_mesh,&
                                        ibool_interfaces_ext_mesh,1,CUDA_AWARE_MPI)


        call wait_sent_request(tab_requests_send_recv_vector_c,ninterface_elastic,NPROC,inum_interfaces_elastic)

       endif !CUDA_AWARE_MPI

    endif  ! Phase_is_inner

  enddo ! Phase_is_inner


 ! multiplies with inverse of mass matrix (note: rmass has been inverted already)
 call kernel_3_a_cuda(Mesh_pointer,deltatover2f,b_deltatover2f)
  
 

end subroutine compute_forces_elastic_GPU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic_GPU(phase_is_inner,num_abs_boundary_faces, &
                        SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                        Mesh_pointer,b_absorb_elastic_bottom_slice,b_absorb_elastic_left_slice,b_absorb_elastic_right_slice,&
                        b_absorb_elastic_top_slice)

  use constants
  use specfem_par, only : nspec_bottom,nspec_left,nspec_top,nspec_right,b_absorb_elastic_left,b_absorb_elastic_right,&
                           b_absorb_elastic_bottom, b_absorb_elastic_top
  implicit none

! communication overlap
  logical :: phase_is_inner

! absorbing boundary surface
  integer :: num_abs_boundary_faces

! adjoint simulations
  integer:: SIMULATION_TYPE
  integer:: NSTEP,it,i,j


  logical:: SAVE_FORWARD


  ! GPU_MODE variables
  integer(kind=8) :: Mesh_pointer
  real(kind=CUSTOM_REAL),dimension(2,NGLLX,nspec_bottom) :: b_absorb_elastic_bottom_slice
  real(kind=CUSTOM_REAL),dimension(2,NGLLX,nspec_left) :: b_absorb_elastic_left_slice
  real(kind=CUSTOM_REAL),dimension(2,NGLLX,nspec_right) :: b_absorb_elastic_right_slice
  real(kind=CUSTOM_REAL),dimension(2,NGLLX,nspec_top) :: b_absorb_elastic_top_slice


  ! checks if anything to do
  if( num_abs_boundary_faces == 0 ) return

if( SIMULATION_TYPE == 3 ) then
    if( phase_is_inner .eqv. .false. ) then

    b_absorb_elastic_bottom_slice(1,:,:)=b_absorb_elastic_bottom(1,:,:,NSTEP-it+1)
    b_absorb_elastic_left_slice(1,:,:)=b_absorb_elastic_left(1,:,:,NSTEP-it+1)
    b_absorb_elastic_right_slice(1,:,:)=b_absorb_elastic_right(1,:,:,NSTEP-it+1)
    b_absorb_elastic_top_slice(1,:,:)=b_absorb_elastic_top(1,:,:,NSTEP-it+1)

    b_absorb_elastic_bottom_slice(2,:,:)=b_absorb_elastic_bottom(3,:,:,NSTEP-it+1)
    b_absorb_elastic_left_slice(2,:,:)=b_absorb_elastic_left(3,:,:,NSTEP-it+1)
    b_absorb_elastic_right_slice(2,:,:)=b_absorb_elastic_right(3,:,:,NSTEP-it+1)
    b_absorb_elastic_top_slice(2,:,:)=b_absorb_elastic_top(3,:,:,NSTEP-it+1)
 
    endif
endif


  call compute_stacey_viscoelastic_cuda(Mesh_pointer,phase_is_inner,b_absorb_elastic_left_slice,&
                   b_absorb_elastic_right_slice,b_absorb_elastic_top_slice,b_absorb_elastic_bottom_slice)

  ! adjoint simulations: stores absorbed wavefield part
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD ) then
    ! writes out absorbing boundary value only when second phase is running
    if( phase_is_inner .eqv. .true. ) then
    b_absorb_elastic_bottom(1,:,:,it) = b_absorb_elastic_bottom_slice(1,:,:)
    b_absorb_elastic_right(1,:,:,it) = b_absorb_elastic_right_slice(1,:,:)
    b_absorb_elastic_top(1,:,:,it) = b_absorb_elastic_top_slice(1,:,:)
    b_absorb_elastic_left(1,:,:,it) = b_absorb_elastic_left_slice(1,:,:)

    b_absorb_elastic_bottom(3,:,:,it) = b_absorb_elastic_bottom_slice(2,:,:)
    b_absorb_elastic_right(3,:,:,it) = b_absorb_elastic_right_slice(2,:,:)
    b_absorb_elastic_top(3,:,:,it) = b_absorb_elastic_top_slice(2,:,:)
    b_absorb_elastic_left(3,:,:,it) = b_absorb_elastic_left_slice(2,:,:)

    endif
  endif


  end subroutine compute_stacey_viscoelastic_GPU

!
!=====================================================================
! for elastic solver on GPU

  subroutine compute_add_sources_viscoelastic_GPU(phase_is_inner, &
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
      call compute_add_sources_el_cuda(Mesh_pointer,phase_is_inner,NSOURCES,it)
  endif ! forward


! adjoint simulations
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! adds adjoint source in this partitions
    if( nadj_rec_local > 0 ) then
      if( it < NSTEP ) then
        call add_sources_el_sim_type_2_or_3(Mesh_pointer,phase_is_inner, NSTEP -it + 1, nadj_rec_local,NSTEP)
      endif ! it
    endif ! nadj_rec_local
  endif !adjoint


! adjoint simulations
  if (SIMULATION_TYPE == 3 .and. nsources_local > 0) then
     call compute_add_sources_el_s3_cuda(Mesh_pointer,phase_is_inner,NSOURCES,NSTEP -it + 1)
  endif ! adjoint


  end subroutine compute_add_sources_viscoelastic_GPU
