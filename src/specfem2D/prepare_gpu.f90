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

  use constants, only: IMAIN,USE_A_STRONG_FORMULATION_FOR_E1
  use specfem_par
  use specfem_par_gpu
  use moving_sources_par, only: init_moving_sources_GPU

  implicit none

  ! local parameters
  real :: free_mb,used_mb,total_mb
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: cosrot_irecf, sinrot_irecf
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmassx,rmassz
  ! PML
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: alphax_store_GPU,alphaz_store_GPU,betax_store_GPU,betaz_store_GPU

  ! checks if anything to do
  if (.not. GPU_MODE) return

  ! GPU_MODE now defined in Par_file
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) "Preparing GPU Fields and Constants on Device."
    call flush_IMAIN()
  endif

  ! safety checks
  if (AXISYM) &
    call stop_the_code('Axisym not implemented on GPU yet.')
  if (NGLLX /= NGLLZ) &
    call stop_the_code('GPU simulations require NGLLX == NGLLZ')

  if (ATTENUATION_VISCOACOUSTIC .and. (.not. USE_A_STRONG_FORMULATION_FOR_E1)) &
    call stop_the_code('GPU simulations require USE_A_STRONG_FORMULATION_FOR_E1 set to true')

  if (ATTENUATION_VISCOELASTIC .and. SIMULATION_TYPE == 3) &
    call stop_the_code('GPU_MODE not supported yet adjoint simulations with attenuation viscoelastic')
  if ((ATTENUATION_VISCOACOUSTIC .or. ATTENUATION_VISCOELASTIC) .and. any_elastic .and. any_acoustic) &
    call stop_the_code('GPU_MODE not supported yet coupled fluid-solid simulations with attenuation')

  if (PML_BOUNDARY_CONDITIONS .and. any_elastic) &
    call stop_the_code('PML on GPU not supported yet for elastic cases')
  if (PML_BOUNDARY_CONDITIONS .and. ATTENUATION_VISCOACOUSTIC) &
    call stop_the_code('PML on GPU not supported yet for viscoacoustic cases')
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
                                SIMULATION_TYPE,P_SV, &
                                nspec_acoustic,nspec_elastic, &
                                ispec_is_acoustic,ispec_is_elastic, &
                                myrank,SAVE_FORWARD, &
                                xir_store_loc, &
                                gammar_store_loc, &
                                NSIGTYPE, seismotypeVec, &
                                nlength_seismogram)

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
                                       ispec_is_anisotropic, &
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
      call prepare_fields_elastic_adj_dev(Mesh_pointer,NDIM*NGLOB_AB,APPROXIMATE_HESS_KL, &
                                          ATTENUATION_VISCOELASTIC,NO_BACKWARD_RECONSTRUCTION)
    endif

    ! frees memory
    deallocate(rmassx,rmassz)
  endif

  if (PML_BOUNDARY_CONDITIONS) then
    call prepare_PML_device(Mesh_pointer, &
                            nspec_PML, &
                            NSPEC_PML_X,NSPEC_PML_Z,NSPEC_PML_XZ, &
                            spec_to_PML_GPU, &
                            abs_normalized, &
                            sngl(ALPHA_MAX_PML), &
                            d0_max, &
                            deltat, &
                            alphax_store_GPU,alphaz_store_GPU, &
                            betax_store_GPU,betaz_store_GPU, &
                            PML_nglob_abs_acoustic,PML_abs_points_acoustic)
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

    if (SIMULATION_TYPE == 3 .and. (.not. NO_BACKWARD_RECONSTRUCTION)) &
      call transfer_b_fields_ac_to_device(NGLOB_AB,b_potential_acoustic, &
                                          b_potential_dot_acoustic,b_potential_dot_dot_acoustic,Mesh_pointer)
  endif

  ! puts elastic initial fields onto GPU
  if (any_elastic) then
    ! transfers forward fields to device with initial values
    call transfer_fields_el_to_device(NDIM*NGLOB_AB,displ_elastic,veloc_elastic,accel_elastic,Mesh_pointer)

    if (SIMULATION_TYPE == 3) then
      ! transfers backward fields to device with initial values
      call transfer_b_fields_to_device(NDIM*NGLOB_AB,b_displ_elastic,b_veloc_elastic,b_accel_elastic,Mesh_pointer)
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

  ! moving sources
  if (SOURCE_IS_MOVING .and. SIMULATION_TYPE == 1) then
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) "  Initializing variables for moving sources ..."
      call flush_IMAIN()
    endif
    call init_moving_sources_GPU()
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
  deallocate(alphax_store_GPU,alphaz_store_GPU,betax_store_GPU,betaz_store_GPU)

  ! synchronizes all processes
  call synchronize_all()

  contains

!----------------------------------------------------------------------

  subroutine init_host_to_dev_variable()

! helper routine for array initialization and time run setup

  use constants, only: IMAIN,IEDGE1,IEDGE2,IEDGE3,IEDGE4,IBOTTOM,IRIGHT,ITOP, &
                       ILEFT,CPML_X_ONLY,CPML_Z_ONLY,CPML_XZ

  implicit none

  ! local parameters
  integer :: i_spec_free,ipoint1D,i,j,ispec,i_source,i_source_local,ispec_PML,ielem
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

  NSPEC_AB = nspec
  NGLOB_AB = nglob

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  number of MPI interfaces                        = ',ninterface
    call flush_IMAIN()
  endif

  ! PML
  if (nspec_PML == 0) PML_BOUNDARY_CONDITIONS = .false.

  if (PML_BOUNDARY_CONDITIONS) then
    ! EB EB : We create spec_to_PML_GPU such that :
    ! spec_to_PML_GPU(ispec) = 0 indicates the element is not in the PML
    ! spec_to_PML_GPU(ispec) \in [1, NSPEC_PML_X] indicates the element is not in the region CPML_X_ONLY
    ! spec_to_PML_GPU(ispec) \in [NSPEC_PML_X + 1, NSPEC_PML_X + NSPEC_PML_Z] indicates the element is in the region CPML_Z_ONLY
    ! spec_to_PML_GPU(ispec) >  NSPEC_PML_X + NSPEC_PML_Z indicates the element is in the region CPML_XZ
    ! Finally, spec_to_PML_GPU(ispec) = ielem, where ielem the local number of the element in the PML
    allocate(spec_to_PML_GPU(nspec))
    spec_to_PML_GPU(:) = 0
    nspec_PML_X = 0
    do ispec = 1,nspec
      if (region_CPML(ispec) == CPML_X_ONLY ) then
        nspec_PML_X = nspec_PML_X+1
        spec_to_PML_GPU(ispec) = nspec_PML_X
      endif
    enddo
    nspec_PML_Z = 0
    do ispec = 1,nspec
      if (region_CPML(ispec) == CPML_Z_ONLY ) then
        nspec_PML_Z = nspec_PML_Z+1
        spec_to_PML_GPU(ispec) = nspec_PML_X + nspec_PML_Z
      endif
    enddo
    nspec_PML_XZ = 0
    do ispec = 1,nspec
      if (region_CPML(ispec) == CPML_XZ ) then
        nspec_PML_XZ = nspec_PML_XZ+1
        spec_to_PML_GPU(ispec) = nspec_PML_X + nspec_PML_Z + nspec_PML_XZ
      endif
    enddo
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  number of PML elements in this process slice    = ',nspec_PML
      write(IMAIN,*) '      elements X only                             = ',nspec_PML_X
      write(IMAIN,*) '      elements Z only                             = ',nspec_PML_Z
      write(IMAIN,*) '      elements XZ                                 = ',nspec_PML_XZ
      call flush_IMAIN()
    endif

    ! Safety check
    if (nspec_PML_X + nspec_PML_Z + nspec_PML_XZ /= nspec_PML) then
      print *,'Error: rank ',myrank,' has invalid number of PML elements ',nspec_PML, &
              ' X ',nspec_PML_X,'Z ',nspec_PML_Z, 'XZ', nspec_PML_XZ
      stop 'Error with the number of PML elements in GPU mode'
    endif

    ! EB EB : We reorganize the arrays abs_normalized and abs_normalized2 that
    ! don't have the correct dimension and new local element numbering
    allocate(abs_normalized_temp(NGLLX,NGLLZ,NSPEC))
    abs_normalized_temp(:,:,:) = abs_normalized(:,:,:)
    deallocate(abs_normalized)
    allocate(abs_normalized(NGLLX,NGLLZ,NSPEC_PML))
    do ispec = 1,nspec
      if (spec_to_PML_GPU(ispec) > 0) abs_normalized(:,:,spec_to_PML_GPU(ispec)) = abs_normalized_temp(:,:,ispec)
    enddo
    deallocate(abs_normalized_temp)

    allocate(alphax_store_GPU(NGLLX,NGLLZ,NSPEC_PML_XZ),alphaz_store_GPU(NGLLX,NGLLZ,NSPEC_PML_XZ), &
             betax_store_GPU(NGLLX,NGLLZ,NSPEC_PML_XZ),betaz_store_GPU(NGLLX,NGLLZ,NSPEC_PML_XZ))
    alphax_store_GPU(:,:,:) = 0.0_CUSTOM_REAL
    alphaz_store_GPU(:,:,:) = 0.0_CUSTOM_REAL
    betax_store_GPU(:,:,:) = 0.0_CUSTOM_REAL
    betaz_store_GPU(:,:,:) = 0.0_CUSTOM_REAL
    do ispec = 1,nspec
      if (region_CPML(ispec) == CPML_XZ) then
        ispec_PML = spec_to_PML(ispec)
        ! element index in range [1,NSPEC_PML_XZ]
        ielem = spec_to_PML_GPU(ispec)-(nspec_PML_X + nspec_PML_Z)
        do j = 1,NGLLZ
          do i = 1,NGLLX
            alphax_store_GPU(i,j,ielem) = sngl(alpha_x_store(i,j,ispec_PML))
            alphaz_store_GPU(i,j,ielem) = sngl(alpha_z_store(i,j,ispec_PML))
            betax_store_GPU(i,j,ielem) = sngl(alpha_x_store(i,j,ispec_PML) + d_x_store(i,j,ispec_PML))
            betaz_store_GPU(i,j,ielem) = sngl(alpha_z_store(i,j,ispec_PML) + d_z_store(i,j,ispec_PML))
           enddo
         enddo
      endif
    enddo
  else
    ! dummy allocations
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
    write(IMAIN,*) '  total number of sources                         = ',NSOURCES
    write(IMAIN,*) '  number of sources in this process slice         = ',nsources_local
    call flush_IMAIN()
  endif

  ! creates arrays to hold only local source entries
  if (nsources_local > 0) then
    allocate(sourcearray_loc(NDIM,NGLLX,NGLLX,nsources_local))
    allocate(source_time_function_loc(nsources_local,NSTEP))
    allocate(ispec_selected_source_loc(nsources_local))
  else
    ! dummy
    allocate(sourcearray_loc(1,1,1,1))
    allocate(source_time_function_loc(1,1))
    allocate(ispec_selected_source_loc(1))
  endif
  source_time_function_loc(:,:) = 0.0_CUSTOM_REAL
  ispec_selected_source_loc(:) = 1
  sourcearray_loc(:,:,:,:) = 0.0_CUSTOM_REAL

  i_source_local = 0
  do i_source = 1, NSOURCES
    if (myrank == islice_selected_source(i_source)) then
      ! source belongs to this process
      i_source_local = i_source_local + 1
      if (i_source_local > nsources_local) stop 'Error with the number of local sources'
      ! stores local sources infos
      source_time_function_loc(i_source_local,:) = source_time_function(i_source,:,1)
      ispec_selected_source_loc(i_source_local)  = ispec_selected_source(i_source)
      if (P_SV) then
        sourcearray_loc(:,:,:,i_source_local) = sourcearrays(:,:,:,i_source)
      else
        sourcearray_loc(1,:,:,i_source_local) = sourcearrays(1,:,:,i_source) ! only single component for SH
      endif
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
!!!!!!! Init to prepare acoustics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! acoustic elements at free surface
  if (myrank == 0) then
    write(IMAIN,*) '  number of acoustic elements at free surface     = ',nelem_acoustic_surface
    call flush_IMAIN()
  endif
  allocate(free_ac_ispec(nelem_acoustic_surface))
  if (nelem_acoustic_surface > 0) then
    ! gets ispec indices for acoustic elements at free surface
    free_ac_ispec(:) = acoustic_surface(1,:)
  endif

  ! free surface
  allocate(free_surface_ij(2,NGLLX,nelem_acoustic_surface))
  free_surface_ij(:,:,:) = 0
  do i_spec_free = 1, nelem_acoustic_surface
    if (acoustic_surface(2,i_spec_free) == acoustic_surface(3,i_spec_free)) then
      do j = 1,NGLLX
        free_surface_ij(1,j,i_spec_free) = acoustic_surface(2,i_spec_free)
      enddo
    else
      j = 1
      do i = acoustic_surface(2,i_spec_free), acoustic_surface(3,i_spec_free)
        free_surface_ij(1,j,i_spec_free) = i
        j = j+1
      enddo
    endif

    if (acoustic_surface(4,i_spec_free) == acoustic_surface(5,i_spec_free)) then
      do j = 1,NGLLX
        free_surface_ij(2,j,i_spec_free) = acoustic_surface(4,i_spec_free)
      enddo
    else
      j = 1
      do i = acoustic_surface(4,i_spec_free), acoustic_surface(5,i_spec_free)
        free_surface_ij(2,j,i_spec_free) = i
        j = j+1
      enddo
    endif
  enddo

  ! coupling surfaces for acoustic-elastic domains
  if (myrank == 0) then
    write(IMAIN,*) '  number of coupled fluid-solid edges             = ',num_fluid_solid_edges
    call flush_IMAIN()
  endif
  allocate(coupling_ac_el_ispec(num_fluid_solid_edges))
  allocate(coupling_ac_el_ij(2,NGLLX,num_fluid_solid_edges))
  allocate(coupling_ac_el_normal(2,NGLLX,num_fluid_solid_edges))
  allocate(coupling_ac_el_jacobian1Dw(NGLLX,num_fluid_solid_edges))
  coupling_ac_el_ispec(:) = 0
  coupling_ac_el_ij(:,:,:) = 0
  coupling_ac_el_normal(:,:,:) = 0.0_CUSTOM_REAL
  coupling_ac_el_jacobian1Dw(:,:) = 0.00_CUSTOM_REAL
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
    write(IMAIN,*)
    write(IMAIN,*) '  init_host_to_dev_variable done successfully'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine init_host_to_dev_variable

  end subroutine prepare_GPU

