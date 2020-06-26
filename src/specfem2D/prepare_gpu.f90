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

  subroutine prepare_GPU()

  use constants, only: IMAIN,APPROXIMATE_HESS_KL,USE_A_STRONG_FORMULATION_FOR_E1,SOURCE_IS_MOVING
  use specfem_par
  use specfem_par_gpu

  implicit none

  ! local parameters
  real :: free_mb,used_mb,total_mb
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: cosrot_irecf, sinrot_irecf
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmassx,rmassz

  ! checks if anything to do
  if (.not. GPU_MODE) return

  ! GPU_MODE now defined in Par_file
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) "Preparing GPU Fields and Constants on Device."
    call flush_IMAIN()
  endif

  ! safety checks
  if (any_elastic .and. (.not. P_SV)) &
    call stop_the_code('Invalid GPU simulation, SH waves not implemented yet. Please use P_SV instead.')
  if (AXISYM) &
    call stop_the_code('Axisym not implemented on GPU yet.')
  if (NGLLX /= NGLLZ) &
    call stop_the_code('GPU simulations require NGLLX == NGLLZ')
  if ( (.not. USE_A_STRONG_FORMULATION_FOR_E1) .and. ATTENUATION_VISCOACOUSTIC) &
    call stop_the_code('GPU simulations require USE_A_STRONG_FORMULATION_FOR_E1 set to true')
  if ( ATTENUATION_VISCOELASTIC .and. SIMULATION_TYPE == 3) &
    call stop_the_code('GPU mode do not support yet adjoint simulations with attenuation viscoelastic')
  if ( (ATTENUATION_VISCOACOUSTIC .or. ATTENUATION_VISCOELASTIC) .and. any_elastic .and. any_acoustic) &
    call stop_the_code('GPU mode do not support yet coupled fluid-solid simulations with attenuation')
  if (PML_BOUNDARY_CONDITIONS .and. any_elastic) &
    call stop_the_code('PML on GPU do not support elastic case yet')
  if (PML_BOUNDARY_CONDITIONS .and. ATTENUATION_VISCOACOUSTIC) &
    call stop_the_code('PML on GPU do not support viscoacoustic case yet')
  if (PML_BOUNDARY_CONDITIONS .and. SIMULATION_TYPE == 3 .and. (.not. NO_BACKWARD_RECONSTRUCTION) ) &
    call stop_the_code('PML on GPU in adjoint mode only work using NO_BACKWARD_RECONSTRUCTION flag')

  ! initializes arrays
  call init_host_to_dev_variable()

  ! check number of purely elastic elements
  if (nspec_elastic /= nspec - nspec_acoustic) then
    print *,'GPU simulation only supported for acoustic and/or elastic domain simulations'
    call stop_the_code('Error GPU simulation')
  endif

  ! Input parameters :

  ! ibool(i,j,ispec)                     : array that converts the index of GLL point (i,j) in elmt ispec
  !                                        from local to global (iglob)
  ! ninterface                           : Number of interfaces that the local partition share with other partitions
  ! max_nibool_interfaces_ext_mesh       : Maximum number of GLL points at an interface
  ! nibool_interfaces_ext_mesh(i)        : Number of GLL points in interface i
  ! ibool_interfaces_ext_mesh(iGGL,i)    : Global index iglob of the ith GLL point (iGLL) of interface i
  ! ispec_is_inner                       : Boolean array. True if the input element is into a partition
  ! sourcearray_loc(i_src,dim,i,j)       : Array of weights containing source intensity at each GLL point (i,j) of the spectral
  !                                        element containing the local source i_src
  ! ispec_selected_source(i)             : Index of the spectral element containing local source i
  ! ispec_selected_rec_loc(i)            : Index of the spectral element containing local receiver i
  ! nrecloc                              : Number of local receivers
  ! nspec_acoustic                       : Number of local acoustic spectral elements

  ! prepares general fields on GPU
  !! JC JC here we will need to add GPU support for the C-PML routines

  call prepare_constants_device(Mesh_pointer, &
                                NGLLX, nspec, nglob, &
                                xix, xiz, gammax, gammaz, &
                                kappastore, mustore, &
                                ibool, &
                                ninterface, max_nibool_interfaces_ext_mesh, &
                                nibool_interfaces_ext_mesh, ibool_interfaces_ext_mesh, &
                                hprime_xx,hprimewgll_xx, &
                                wxgll, &
                                STACEY_ABSORBING_CONDITIONS, &
                                PML_BOUNDARY_CONDITIONS, &
                                ispec_is_inner, &
                                nsources_local, &
                                sourcearray_loc,source_time_function_loc, &
                                NSTEP,ispec_selected_source_loc, &
                                ispec_selected_rec_loc, &
                                nrecloc, &
                                cosrot_irecf,sinrot_irecf, &
                                SIMULATION_TYPE, &
                                nspec_acoustic,nspec_elastic, &
                                myrank,SAVE_FORWARD, &
                                xir_store_loc, &
                                gammar_store_loc, &
                                NSIGTYPE, seismotypeVec, &
                                NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos)

  ! Input parameters :

  ! rmass_inverse_acoustic                 : Inverse of the acoustic mass matrix (size nglob)
  !                                          (nglob_acoustic = nglob if it exists acoustic elements)
  ! num_phase_ispec_acoustic               : Max between the number of inner acoustic spectral elements and outer
  ! phase_ispec_inner_acoustic(i,j)        : ith acoustic spectral element. Inner if j=2 outer if j=1
  ! acoustic(i)                            : True if spectral element i is acoustic
  ! nelem_acoustic_surface                 : Number of spectral elements on an acoustic free surface
  ! free_ac_ispec                          : Index of ith acoustic spectral element on free surfaces
  ! free_surface_ij(i,j,ispec)             : ith coordinate of jth GLL point of the spectral element ispec located on a free surface
  ! b_reclen_potential                     : Size in bytes taken by b_nelem_acoustic_surface * GLLX
  ! any_elastic                            : True if it exists elastic elements in the mesh
  ! num_fluid_solid_edges                  : Number of spectral elements on the elasto-acoustic border
  ! coupling_ac_el_ispec                   : Array containing spectral element indices on the elasto-acoustic border
  ! coupling_ac_el_ij                      : Local coordinates of GLL points on the elasto-acoustic border
  ! coupling_ac_el_normal(i,j,ispec)       : ith coordinate of normal vector at GLL point j in border element ispec
  ! coupling_ac_el_jacobian1Dw(i,ispec)    : Weighted jacobian matrix at ith GLL point in border element ispec

  ! prepares fields on GPU for acoustic simulations
  if (any_acoustic) then
    call prepare_fields_acoustic_device(Mesh_pointer, &
                                        rmass_inverse_acoustic,rhostore,kappastore, &
                                        num_phase_ispec_acoustic,phase_ispec_inner_acoustic, &
                                        ispec_is_acoustic, &
                                        nelem_acoustic_surface, &
                                        free_ac_ispec,free_surface_ij, &
                                        any_elastic, num_fluid_solid_edges, &
                                        coupling_ac_el_ispec,coupling_ac_el_ij, &
                                        coupling_ac_el_normal,coupling_ac_el_jacobian1Dw, &
                                        ninterface_acoustic,inum_interfaces_acoustic,ATTENUATION_VISCOACOUSTIC, &
                                        A_newmark_e1_sf,B_newmark_e1_sf,NO_BACKWARD_RECONSTRUCTION,no_backward_acoustic_buffer)

    if (SIMULATION_TYPE == 3) then
      ! safety check
      if (APPROXIMATE_HESS_KL) then
        call stop_the_code('Sorry, approximate acoustic Hessian kernels not yet fully implemented for GPU simulations!')
      endif
      call prepare_fields_acoustic_adj_dev(Mesh_pointer,APPROXIMATE_HESS_KL, &
                                           ATTENUATION_VISCOACOUSTIC,NO_BACKWARD_RECONSTRUCTION)
    endif
  endif

  ! Input parameters :

  ! rmass_inverse_elastic          : Inverse elastic mass matrix (size nglob_acoustic)
  !                                  (nglob_acoustic = nglob if there are acoustic elements)
  ! num_phase_ispec_elastic        : Max between the number of inner elastic spectral elements and outer
  ! phase_ispec_inner_elastic(i,j) : ith elastic spectral element. Inner if j=2 outer if j=1
  ! elastic(i)                     : True if spectral element i is elastic

  ! prepares fields on GPU for elastic simulations
  if (any_elastic) then
    ! temporary mass matrices
    allocate(rmassx(nglob_elastic),rmassz(nglob_elastic))
    rmassx(:) = rmass_inverse_elastic(1,:)
    rmassz(:) = rmass_inverse_elastic(2,:)

    call prepare_fields_elastic_device(Mesh_pointer, &
                                       rmassx,rmassz, &
                                       num_phase_ispec_elastic,phase_ispec_inner_elastic, &
                                       ispec_is_elastic, &
                                       any_anisotropy, &
                                       c11store,c12store,c13store, &
                                       c15store,c23store, &
                                       c25store,c33store,c35store,c55store, &
                                       ninterface_elastic,inum_interfaces_elastic, &
                                       ATTENUATION_VISCOELASTIC, &
                                       A_newmark_nu2,B_newmark_nu2,A_newmark_nu1,B_newmark_nu1)


    if (SIMULATION_TYPE == 3) then
      ! safety check
      if (APPROXIMATE_HESS_KL) then
        call stop_the_code('Sorry, approximate elastic Hessian kernels not yet fully implemented for GPU simulations!')
      endif
      call prepare_fields_elastic_adj_dev(Mesh_pointer,NDIM*NGLOB_AB,APPROXIMATE_HESS_KL)
    endif

    ! frees memory
    deallocate(rmassx,rmassz)
  endif

  if (PML_BOUNDARY_CONDITIONS) then
    call prepare_PML_device(Mesh_pointer, &
                            nspec_PML, &
                            NSPEC_PML_X, &
                            NSPEC_PML_Z, &
                            NSPEC_PML_XZ, &
                            spec_to_PML_GPU, &
                            abs_normalized, &
                            sngl(ALPHA_MAX_PML), &
                            d0_max, &
                            sngl(deltat), &
                            alphax_store_GPU, &
                            alphaz_store_GPU, &
                            betax_store_GPU, &
                            betaz_store_GPU)
  endif

! abs_boundary_ispec                     : Array containing spectral element indices in absorbing areas
! abs_boundary_ij(i,j,ispecabs)          : ith local coordinate of jth GLL point of the absorbing element ispecabs
! abs_boundary_normal(i,j,ispecabs)      : ith coordinate of normal vector at GLL point j in absorbing element ispecabs
! abs_boundary_jacobian1Dw(i,ispecabs)   : Weighted jacobian matrix at ith GLL point in absorbing element jspecabs
! num_abs_boundary_faces                 : Number of absorbing elements
! edge_abs(ispecabs)                     : Index of edge (1=bottom, 2=right, 3=top, 4=left) corresponding to absorbing elmt ispecabs
! ib_left                                : Correpondance between the global absorbing element index and its index on the edge

  if (STACEY_ABSORBING_CONDITIONS) then
    call prepare_Stacey_device(Mesh_pointer, &
                               any_acoustic,any_elastic, &
                               rho_vpstore,rho_vsstore, &
                               nspec_bottom, &
                               nspec_left, &
                               nspec_right, &
                               nspec_top, &
                               abs_boundary_ispec, abs_boundary_ij, &
                               abs_boundary_normal, &
                               abs_boundary_jacobian1Dw, &
                               num_abs_boundary_faces, &
                               edge_abs, &
                               ib_bottom, &
                               ib_left, &
                               ib_right, &
                               ib_top)
  endif

  ! prepares fields on GPU for poroelastic simulations
  if (any_poroelastic) then
    call stop_the_code('todo poroelastic simulations on GPU')
  endif

  ! prepares needed receiver array for adjoint runs
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) &
    call prepare_sim2_or_3_const_device(Mesh_pointer,nrecloc,source_adjoint,NSTEP)

  ! synchronizes processes
  call synchronize_all()

  ! puts acoustic initial fields onto GPU
  if (any_acoustic) then
    call transfer_fields_ac_to_device(NGLOB_AB,potential_acoustic, &
                                      potential_dot_acoustic,potential_dot_dot_acoustic,Mesh_pointer)

    if (SIMULATION_TYPE == 3 .and. .not. NO_BACKWARD_RECONSTRUCTION) &
      call transfer_b_fields_ac_to_device(NGLOB_AB,b_potential_acoustic, &
                                          b_potential_dot_acoustic,b_potential_dot_dot_acoustic,Mesh_pointer)
  endif

  ! puts elastic initial fields onto GPU
  if (any_elastic) then
    ! prepares wavefields for transfering
    allocate(tmp_displ_2D(NDIM,nglob_elastic), &
             tmp_veloc_2D(NDIM,nglob_elastic), &
             tmp_accel_2D(NDIM,nglob_elastic))

    tmp_displ_2D(:,:) = displ_elastic(:,:)
    tmp_veloc_2D(:,:) = veloc_elastic(:,:)
    tmp_accel_2D(:,:) = accel_elastic(:,:)

    ! transfers forward fields to device with initial values
    call transfer_fields_el_to_device(NDIM*NGLOB_AB,tmp_displ_2D,tmp_veloc_2D,tmp_accel_2D,Mesh_pointer)

    if (SIMULATION_TYPE == 3) then
      tmp_displ_2D(:,:) = b_displ_elastic(:,:)
      tmp_veloc_2D(:,:) = b_veloc_elastic(:,:)
      tmp_accel_2D(:,:) = b_accel_elastic(:,:)

      ! transfers backward fields to device with initial values
      call transfer_b_fields_to_device(NDIM*NGLOB_AB,tmp_displ_2D,tmp_veloc_2D,tmp_accel_2D,Mesh_pointer)
    endif
  endif

  ! allocates arrays for MPI transfers
  allocate(request_send_recv_scalar_gpu(2*ninterface))
  allocate(b_request_send_recv_scalar_gpu(2*ninterface))

  allocate(request_send_recv_vector_gpu(2*ninterface))
  allocate(b_request_send_recv_vector_gpu(2*ninterface))

  allocate(buffer_send_scalar_gpu(max_nibool_interfaces_ext_mesh,ninterface))
  allocate(b_buffer_send_scalar_gpu(max_nibool_interfaces_ext_mesh,ninterface))

  allocate(buffer_recv_scalar_gpu(max_nibool_interfaces_ext_mesh,ninterface))
  allocate(b_buffer_recv_scalar_gpu(max_nibool_interfaces_ext_mesh,ninterface))

  allocate(buffer_send_vector_gpu(2,max_nibool_interfaces_ext_mesh,ninterface))
  allocate(b_buffer_send_vector_gpu(2,max_nibool_interfaces_ext_mesh,ninterface))

  allocate(buffer_recv_vector_gpu(2,max_nibool_interfaces_ext_mesh,ninterface))
  allocate(b_buffer_recv_vector_gpu(2,max_nibool_interfaces_ext_mesh,ninterface))

  if (SOURCE_IS_MOVING .and. SIMULATION_TYPE == 1) then
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) "  Initializing variables for moving sources ..."
      call flush_IMAIN()
    endif
    call init_moving_source()
    if (myrank == 0) then
      write(IMAIN,*) "  Done"
      call flush_IMAIN()
    endif
  endif
  ! synchronizes processes
  call synchronize_all()

  ! outputs GPU usage to files for all processes
  call output_free_device_memory(myrank)

  ! outputs usage for main process
  if (myrank == 0) then
    write(IMAIN,*) "Initialization done successfully"
    write(IMAIN,*)
    call get_free_device_memory(free_mb,used_mb,total_mb)

    write(IMAIN,*)
    write(IMAIN,*) "GPU usage: free  =",free_mb," MB",nint(free_mb/total_mb*100.0),"%"
    write(IMAIN,*) "           used  =",used_mb," MB",nint(used_mb/total_mb*100.0),"%"
    write(IMAIN,*) "           total =",total_mb," MB",nint(total_mb/total_mb*100.0),"%"
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! frees temporary arrays
  deallocate(cosrot_irecf,sinrot_irecf)

  ! synchronizes all processes
  call synchronize_all()

  contains

!----------------------------------------------------------------------

  subroutine init_host_to_dev_variable()

! helper routine for array initialization and time run setup

  use constants, only: IMAIN,IEDGE1,IEDGE2,IEDGE3,IEDGE4,IBOTTOM,IRIGHT,ITOP,ILEFT,CPML_X_ONLY,CPML_Z_ONLY,CPML_XZ

  implicit none

  ! local parameters
  integer :: i_spec_free, ipoint1D, i, j, k, ispec, i_source
  integer :: ispec_acoustic,ispec_elastic,iedge_acoustic,iedge_elastic
  integer :: inum
  real(kind=CUSTOM_REAL) :: zxi,xgamma,jacobian1D
  real(kind=CUSTOM_REAL) :: xxi,zgamma
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: abs_normalized_temp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Initialize variables for subroutine prepare_constants_device
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  initializing arrays for GPU:'
    call flush_IMAIN()
  endif

  ! converts to CUSTOM_REAL
  deltatf = sngl(deltat)
  deltatover2f = sngl(deltatover2)
  deltatsquareover2f = sngl(deltatsquareover2)

  b_deltatf = sngl(b_deltat)
  b_deltatover2f = sngl(b_deltatover2)
  b_deltatsquareover2f = sngl(b_deltatsquareover2)

  NSPEC_AB = nspec
  NGLOB_AB = nglob

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  number of MPI interfaces                        = ',ninterface
    write(IMAIN,*) '  number of acoustic elements at free surface     = ',nelem_acoustic_surface
    call flush_IMAIN()
  endif

  ! acoustic elements at free surface
  allocate(free_ac_ispec(nelem_acoustic_surface))
  if (nelem_acoustic_surface > 0) then
    ! gets ispec indices for acoustic elements at free surface
    free_ac_ispec(:) = acoustic_surface(1,:)
  endif

  if (PML_BOUNDARY_CONDITIONS) then
    ! EB EB : We create spec_to_PML_GPU such that :
    ! spec_to_PML_GPU(ispec) = 0 indicates the element is not in the PML
    ! spec_to_PML_GPU(ispec) \in [1, NSPEC_PML_X] indicates the element is not in the region CPML_X_ONLY
    ! spec_to_PML_GPU(ispec) \in [NSPEC_PML_X + 1, NSPEC_PML_X + NSPEC_PML_Z] indicates the element is in the region CPML_Z_ONLY
    ! spec_to_PML_GPU(ispec) >  NSPEC_PML_X + NSPEC_PML_Z indicates the element is in the region CPML_XZ
    ! Finally, spec_to_PML_GPU(ispec) = ispec_pml, where ispec_pml the local
    ! number of the element in the PML
    allocate(spec_to_PML_GPU(nspec))
    spec_to_PML_GPU(:) = 0
    nspec_PML_X = 0
    do ispec= 1,nspec
      if (region_CPML(ispec) == CPML_X_ONLY ) then
        nspec_PML_X = nspec_PML_X+1
        spec_to_PML_GPU(ispec) = nspec_PML_X
      endif
    enddo
    nspec_PML_Z = 0
    do ispec= 1,nspec
      if (region_CPML(ispec) == CPML_Z_ONLY ) then
        nspec_PML_Z = nspec_PML_Z+1
        spec_to_PML_GPU(ispec) = nspec_PML_X + nspec_PML_Z
      endif
    enddo
    nspec_PML_XZ = 0
    do ispec= 1,nspec
      if (region_CPML(ispec) == CPML_XZ ) then
        nspec_PML_XZ = nspec_PML_XZ+1
        spec_to_PML_GPU(ispec) = nspec_PML_X + nspec_PML_Z + nspec_PML_XZ
      endif
    enddo
    ! Safety check
    if (nspec_PML_X + nspec_PML_Z + nspec_PML_XZ /= nspec_PML) &
      stop 'Error with the number of PML elements in GPU mode'

    ! EB EB : We reorganize the arrays abs_normalized and abs_normalized2 that
    ! don't have the correct dimension and new local element numbering
    allocate(abs_normalized_temp(NGLLX,NGLLZ,NSPEC))
    abs_normalized_temp = abs_normalized
    deallocate(abs_normalized)
    allocate(abs_normalized(NGLLX,NGLLZ,NSPEC_PML))
    do ispec= 1,nspec
     if (spec_to_PML_GPU(ispec) > 0) abs_normalized(:,:,spec_to_PML_GPU(ispec)) = abs_normalized_temp(:,:,ispec)
    enddo
    deallocate(abs_normalized_temp)

    allocate(alphax_store_GPU(NGLLX,NGLLZ,NSPEC_PML_XZ),alphaz_store_GPU(NGLLX,NGLLZ,NSPEC_PML_XZ), &
             betax_store_GPU(NGLLX,NGLLZ,NSPEC_PML_XZ),betaz_store_GPU(NGLLX,NGLLZ,NSPEC_PML_XZ))
    do ispec= 1,nspec
      if (region_CPML(ispec) == CPML_XZ) then
        do j=1,NGLLZ
          do i = 1,NGLLX
            alphax_store_GPU(i,j,spec_to_PML_GPU(ispec)-(nspec_PML_X + nspec_PML_Z)) = sngl(alpha_x_store(i,j,spec_to_PML(ispec)))
            alphaz_store_GPU(i,j,spec_to_PML_GPU(ispec)-(nspec_PML_X + nspec_PML_Z)) = sngl(alpha_z_store(i,j,spec_to_PML(ispec)))
            betax_store_GPU(i,j,spec_to_PML_GPU(ispec)-(nspec_PML_X + nspec_PML_Z)) = sngl(alpha_x_store(i,j,spec_to_PML(ispec)) + &
              d_x_store(i,j,spec_to_PML(ispec)))
            betaz_store_GPU(i,j,spec_to_PML_GPU(ispec)-(nspec_PML_X + nspec_PML_Z)) = sngl(alpha_z_store(i,j,spec_to_PML(ispec)) + &
              d_z_store(i,j,spec_to_PML(ispec)))
           enddo
         enddo
      endif
    enddo
  else
    allocate(spec_to_PML_GPU(1))
    allocate(alphax_store_GPU(1,1,1),alphaz_store_GPU(1,1,1), &
             betax_store_GPU(1,1,1),betaz_store_GPU(1,1,1))
  endif ! PML_BOUNDARY_CONDITIONS

  ! sources
  ! counts sources in this process slice
  nsources_local = 0
  do i = 1, NSOURCES
    if (myrank == islice_selected_source(i)) then
      nsources_local = nsources_local + 1
    endif
  enddo

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  number of sources                               = ',NSOURCES
    write(IMAIN,*) '  number of sources in this process slice         = ',nsources_local
    call flush_IMAIN()
  endif

  allocate(source_time_function_loc(nsources_local,NSTEP))
  allocate(ispec_selected_source_loc(nsources_local))

  j = 0
  do i = 1, NSOURCES
    if (myrank == islice_selected_source(i)) then
      if (j > nsources_local) stop 'Error with the number of local sources'
      j = j + 1
      source_time_function_loc(j,:) = source_time_function(i,:,1)
      ispec_selected_source_loc(j)  = ispec_selected_source(i)
    endif
  enddo

  if (nsources_local > 0) then
    allocate(sourcearray_loc(nsources_local,NDIM,NGLLX,NGLLX))
  else
    allocate(sourcearray_loc(1,1,1,1))
  endif

  k = 0
  do i_source = 1,NSOURCES
    if (myrank == islice_selected_source(i_source)) then
      ! source belongs to this process
      k = k + 1
      sourcearray_loc(k,:,:,:) = sourcearrays(i_source,:,:,:)
    endif
  enddo

  ! receivers
  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  number of local receivers in this process slice = ',nrecloc
    call flush_IMAIN()
  endif

  ! converts to CUSTOM_REAL arrays
  allocate(cosrot_irecf(nrecloc), &
           sinrot_irecf(nrecloc))
  do i = 1,nrecloc
    cosrot_irecf(i) = sngl(cosrot_irec(i))
    sinrot_irecf(i) = sngl(sinrot_irec(i))
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Init pour prepare acoustique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! free surface
  allocate(free_surface_ij(2,NGLLX,nelem_acoustic_surface))

  do i_spec_free = 1, nelem_acoustic_surface
    if (acoustic_surface(2,i_spec_free) == acoustic_surface(3,i_spec_free)) then
      do j =1,NGLLX
        free_surface_ij(1,j,i_spec_free) = acoustic_surface(2,i_spec_free)
      enddo
    else
      j=1
      do i = acoustic_surface(2,i_spec_free), acoustic_surface(3,i_spec_free)
        free_surface_ij(1,j,i_spec_free) = i
        j=j+1
      enddo
    endif

    if (acoustic_surface(4,i_spec_free) == acoustic_surface(5,i_spec_free)) then
      do j =1,NGLLX
        free_surface_ij(2,j,i_spec_free) = acoustic_surface(4,i_spec_free)
      enddo
    else
      j=1
      do i = acoustic_surface(4,i_spec_free), acoustic_surface(5,i_spec_free)
        free_surface_ij(2,j,i_spec_free) = i
        j=j+1
      enddo
    endif
  enddo

  ! coupling surfaces for acoustic-elastic domains
  allocate(coupling_ac_el_ispec(num_fluid_solid_edges))
  allocate(coupling_ac_el_ij(2,NGLLX,num_fluid_solid_edges))
  allocate(coupling_ac_el_normal(2,NGLLX,num_fluid_solid_edges))
  allocate(coupling_ac_el_jacobian1Dw(NGLLX,num_fluid_solid_edges))

  do inum = 1,num_fluid_solid_edges
    ! get the edge of the acoustic element
    ispec_acoustic = fluid_solid_acoustic_ispec(inum)
    iedge_acoustic = fluid_solid_acoustic_iedge(inum)
    coupling_ac_el_ispec(inum) = ispec_acoustic

    ! get the corresponding edge of the elastic element
    ispec_elastic = fluid_solid_elastic_ispec(inum)
    iedge_elastic = fluid_solid_elastic_iedge(inum)

    ! implement 1D coupling along the edge
    do ipoint1D = 1,NGLLX

      ! get point values for the elastic side, which matches our side in the inverse direction
      coupling_ac_el_ij(1,ipoint1D,inum) = ivalue(ipoint1D,iedge_acoustic)
      coupling_ac_el_ij(2,ipoint1D,inum) = jvalue(ipoint1D,iedge_acoustic)

      i = ivalue(ipoint1D,iedge_acoustic)
      j = jvalue(ipoint1D,iedge_acoustic)

      if (iedge_acoustic == ITOP) then
        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        coupling_ac_el_normal(1,ipoint1D,inum) = - zxi / jacobian1D
        coupling_ac_el_normal(2,ipoint1D,inum) = + xxi / jacobian1D
        coupling_ac_el_jacobian1Dw(ipoint1D,inum) = jacobian1D * wxgll(i)

      else if (iedge_acoustic == IBOTTOM) then
        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        coupling_ac_el_normal(1,ipoint1D,inum) = + zxi / jacobian1D
        coupling_ac_el_normal(2,ipoint1D,inum) = - xxi / jacobian1D
        coupling_ac_el_jacobian1Dw(ipoint1D,inum) = jacobian1D * wxgll(i)

      else if (iedge_acoustic == ILEFT) then
        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        coupling_ac_el_normal(1,ipoint1D,inum) = - zgamma / jacobian1D
        coupling_ac_el_normal(2,ipoint1D,inum) = + xgamma / jacobian1D
        coupling_ac_el_jacobian1Dw(ipoint1D,inum) = jacobian1D * wzgll(j)

      else if (iedge_acoustic == IRIGHT) then
        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        coupling_ac_el_normal(1,ipoint1D,inum) = + zgamma / jacobian1D
        coupling_ac_el_normal(2,ipoint1D,inum) = - xgamma / jacobian1D
        coupling_ac_el_jacobian1Dw(ipoint1D,inum) = jacobian1D * wzgll(j)
      endif
    enddo
  enddo

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  init_host_to_dev_variable done successfully'
    call flush_IMAIN()
  endif

  end subroutine init_host_to_dev_variable

  end subroutine prepare_GPU


!
!----------------------------------------------------------------------
!

! note: putting this subroutine as a contains-subroutine into prepare_GPU(),
!       gnu compilers (>version 9) will complain about unused dummy arguments:
!         Warning: Unused dummy argument '_formal_0' at (1) [-Wunused-dummy-argument]
!       in every subroutine call which includes some arguments.
!
!       we thus move this routine out of prepare_GPU() and have it stand alone as subroutine.
!       in this case, the gnu compilers will be happy again (until somebody fixes the gnu compilers...)

  subroutine init_moving_source()

  ! Compute the sourcearrays for all timesteps and send them to device

  use constants, only: NGLLX,NGLLZ,TINYVAL,IMAIN

  use specfem_par
  use specfem_par_gpu

  implicit none

  ! Local variables
  logical :: writeMovingDatabases, file_exists
  character(len=MAX_STRING_LEN) :: pathToMovingDatabase, outputname
  integer :: i_source,i,j,k,ispec,it_l,i_stage_loc,ier
  double precision :: hlagrange
  double precision :: xminSource,vSource,time_val,t_used
  ! Single source array
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: sourcearray_temp

  allocate(ispec_selected_source_loc_moving(1,1))
  allocate(sourcearray_loc_moving(1,1,1,1,1))
  nsources_local_moving(:) = 0
  sourcearrays_moving(:,:,:,:,:) = 0
  file_exists = .false.

  writeMovingDatabases = .true. ! TODO Turn on this if the generation of moving source databases takes a long time
  ! Then you can do it in two steps, first you run the code and it write the databases to file. Then you rerun the code
  ! and it will read the databases from there directly saving a lot of time.

  if (writeMovingDatabases) then
    ! saves data in a binary file
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_data.bin'
    pathToMovingDatabase = '/path/to/wherever/directory/you/want/'  ! Must end with '/'
    inquire(file = trim(pathToMovingDatabase)//outputname, exist=file_exists)
  endif

  if ((.not. file_exists) .or. (.not. writeMovingDatabases)) then  ! Generate the moving source databases (can be expensive)
    if (myrank == 0) then
      write(IMAIN,*) '    Generating moving source databases ...'
      write(IMAIN,*) '      (if this step takes a long time in your case you may want to turn writeMovingDatabases to .true. .'
      write(IMAIN,*) '      In this case goto init_moving_source. Read the comments also)'
      call flush_IMAIN()
    endif
    do it_l = 1, NSTEP  ! Loop over the time steps to precompute all source positions
      do i_stage_loc = 1, stage_time_scheme  ! For now stage_time_scheme = 1 only because only Newmark time scheme are supported

        xminSource = -200.0d0 ! m
        vSource = 200.0d0 ! m/s

        if (time_stepping_scheme == 1) then  ! Only Newmark time scheme is supported
          ! Newmark
          time_val = (it_l-1)*deltat
        else
          call exit_MPI(myrank,'Not implemented! (546999)')
        endif
        ! moves and re-locates sources along x-axis
        do i_source = 1,NSOURCES
          if (abs(source_time_function(i_source,it_l,i_stage_loc)) > TINYVAL) then
            t_used = (time_val-t0-tshift_src(i_source))

            x_source(i_source) = xminSource + vSource*t_used !time_val?

            ! collocated force source (most expensive step)
            call locate_source(ibool,coord,nspec,nglob,xigll,zigll, &
                               x_source(i_source),z_source(i_source), &
                               ispec_selected_source_moving(i_source, it_l),islice_selected_source(i_source), &
                               NPROC,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,ngnod,npgeo, &
                               iglob_source(i_source),.true.)

            call lagrange_any(xi_source(i_source),NGLLX,xigll,hxis,hpxis)
            call lagrange_any(gamma_source(i_source),NGLLZ,zigll,hgammas,hpgammas)

            ! stores Lagrangians for source
            hxis_store(i_source,:) = hxis(:)
            hgammas_store(i_source,:) = hgammas(:)

            sourcearray_temp(:,:,:) = 0._CUSTOM_REAL

            ! element containing source
            ispec = ispec_selected_source_moving(i_source, it_l)
            if (myrank == islice_selected_source(i_source)) then
              ! computes source arrays
              if (source_type(i_source) == 1) then
                ! collocated force source
                do j = 1,NGLLZ
                  do i = 1,NGLLX
                    hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
                    ! source element is acoustic
                    if (ispec_is_acoustic(ispec)) then
                      sourcearray_temp(:,i,j) = hlagrange
                    else
                      call exit_MPI(myrank,'Moving source not in acoustic element (GPU)')
                    endif
                  enddo
                enddo
              else
                call exit_MPI(myrank,'Moving source with source_type != 1 (GPU)')
              endif

              ! stores sourcearray for all sources
              sourcearrays_moving(i_source,:,:,:,it_l) = sourcearray_temp(:,:,:)

            endif
          endif
        enddo

        ! counts sources in this process slice for this time step
        nsources_local = 0
        do i = 1, NSOURCES
          if (myrank == islice_selected_source(i)) then
            nsources_local = nsources_local + 1
          endif
        enddo
        nsources_local_moving(it_l) = nsources_local
      enddo  ! stage_time_scheme
    enddo  ! NSTEP

    deallocate(ispec_selected_source_loc_moving)
    max_nsources_local_moving = maxval(nsources_local_moving)
    allocate(ispec_selected_source_loc_moving(max_nsources_local_moving, NSTEP))
    ispec_selected_source_loc_moving(:,:) = 0

    do it_l = 1, NSTEP
      j = 0
      do i = 1, NSOURCES
        if (myrank == islice_selected_source(i)) then
          if (j > nsources_local_moving(it_l)) call stop_the_code('Error with the number of local sources')
          j = j + 1
          ispec_selected_source_loc_moving(j, it_l)  = ispec_selected_source_moving(i, it_l)
        endif
      enddo
    enddo

    deallocate(sourcearray_loc_moving)
    if (max_nsources_local_moving > 0) then
      allocate(sourcearray_loc_moving(max_nsources_local_moving,NDIM,NGLLX,NGLLX,NSTEP))
    else
      allocate(sourcearray_loc_moving(1,1,1,1,1))
    endif
    do it_l = 1, NSTEP
      k = 0
      do i_source = 1,NSOURCES
        if (myrank == islice_selected_source(i_source)) then
          ! source belongs to this process
          k = k + 1
          sourcearray_loc_moving(k,:,:,:,it_l) = sourcearrays_moving(i_source,:,:,:,it_l)
        endif
      enddo
    enddo

    if (writeMovingDatabases) then  ! Write this to file
      ! user output
      if (myrank == 0) then
        write(IMAIN,*) '  Moving source databases will be written in',trim(pathToMovingDatabase)//trim(outputname),' ...'
        call flush_IMAIN()
      endif
      open(unit = 1698, file = trim(pathToMovingDatabase)//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
      if (ier /= 0) call stop_the_code('  Error writing moving source data file to disk !')
      write(1698) max_nsources_local_moving
      write(1698) NSTEP
      write(1698) nsources_local_moving
      write(1698) sourcearray_loc_moving
      write(1698) ispec_selected_source_loc_moving
      ! user output
      close(1698)
      write(IMAIN,*) '  Moving source databases have been written in ', trim(pathToMovingDatabase)//trim(outputname)
      call synchronize_all()
      if (myrank == 0) then
        write(IMAIN,*) '  Just rerun the code now, they will be read there'
        call flush_IMAIN()
      endif
      call exit_MPI(myrank, '  Terminating ...')
    endif
  else  ! Read the file
    ! Read setup data from a binary file
    write(IMAIN,*) '  Reading moving source databases from file: ', trim(pathToMovingDatabase)//trim(outputname)
    open(unit = 1698,file = trim(pathToMovingDatabase)//outputname,status='old',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_data.bin')

    read(1698) max_nsources_local_moving
    read(1698) NSTEP
    read(1698) nsources_local_moving
    deallocate(ispec_selected_source_loc_moving)
    deallocate(sourcearray_loc_moving)
    allocate(ispec_selected_source_loc_moving(max_nsources_local_moving, NSTEP))
    allocate(sourcearray_loc_moving(max_nsources_local_moving,NDIM,NGLLX,NGLLX,NSTEP))

    read(1698) sourcearray_loc_moving
    read(1698) ispec_selected_source_loc_moving
    close(1698)
  endif

  ! Send variables to device:
  call prepare_moving_source_cuda(Mesh_pointer, nsources_local_moving, max_nsources_local_moving, &
                                  sourcearray_loc_moving, ispec_selected_source_loc_moving, NSTEP)

  end subroutine init_moving_source
