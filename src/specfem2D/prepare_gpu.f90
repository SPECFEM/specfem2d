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

  use constants, only: IMAIN,APPROXIMATE_HESS_KL,USE_A_STRONG_FORMULATION_FOR_E1
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
  if (any_elastic .and. (.not. P_SV)) call stop_the_code( &
'Invalid GPU simulation, SH waves not implemented yet. Please use P_SV instead.')
  if (PML_BOUNDARY_CONDITIONS ) call stop_the_code('PML not implemented on GPU mode. Please use Stacey instead.')
  if (AXISYM) call stop_the_code('Axisym not implemented on GPU yet.')
  if (NGLLX /= NGLLZ) call stop_the_code('GPU simulations require NGLLX == NGLLZ')
  if ( (.not. USE_A_STRONG_FORMULATION_FOR_E1) .and. ATTENUATION_VISCOACOUSTIC) call stop_the_code( &
    'GPU simulations require USE_A_STRONG_FORMULATION_FOR_E1 set to true')
  if ( ATTENUATION_VISCOELASTIC .and. SIMULATION_TYPE == 3) call stop_the_code( &
    'GPU mode do not support yet adjoint simulations with attenuation viscoelastic')
  if ( (ATTENUATION_VISCOACOUSTIC .or. ATTENUATION_VISCOELASTIC) .and. any_elastic .and. any_acoustic) call stop_the_code( &
    'GPU mode do not support yet coupled fluid-solid simulations with attenuation')
  ! initializes arrays
  call init_host_to_dev_variable()

  ! check number of purely elastic elements
  if (nspec_elastic /= nspec - nspec_acoustic) then
    print *,'GPU simulation only supported for acoustic and/or elastic domain simulations'
    call stop_the_code('Error GPU simulation')
  endif

!!!!!!!!!!! Parametres fournis

! ibool(i,j,ispec)                       : convertisseur numero du point GLL local (i,j) de l'element ispec => global (iglob)
! ninterface                             : nombre d'interfaces de la partition locale avec les autres partitions
! max_nibool_interfaces_ext_mesh         : nombre maximum de points GLL contenus sur une interface
! nibool_interfaces_ext_mesh(i)          : nombre de points GLL contenus sur l'interface i
! ibool_interfaces_ext_mesh(iGGL,i)      : numero global iglob du ieme point GLL (iGLL) de l'interface i
! numabs                                 : tableau des elements spectraux situes en zone absorbante
! abs_boundary_ij(i,j,ispecabs)          : coordonnee locale i de j eme point GLL de l'element absorbant ispecabs
! abs_boundary_normal(i,j,ispecabs)      : i eme coordonnee du vecteur normal du j eme point GLL de l'element absorbant ispecabs
! abs_boundary_jacobian1Dw(i,ispecabs)   : i eme jacobienne ponderee de l'element absorbant jspecabs
! nelemabs                               : nombre d'elements absorbant
! cote_abs(ispecabs)                     : numero du cote (1=b, 2=r, 3=t, 4=l ) auquel appartient l'element absorbant ispecabs
! ib_left                                : correspondance entre le numero d'element absorbant global et son numero sur le cote
! ispec_is_inner                         : booleen vrai si l'element est a l'interieur d'une partition
! sourcearray_loc(i_src,dim,i,j)         : tableau de ponderation de l'intensite de la source pour chaque point GLL (i,j)
!                                          de l'element spectral qui contient la source locale i_src
! ispec_selected_source(i)               : numero d'element spectral contenant la source locale i
! ispec_selected_rec_loc(i)              : numero d'element spectral du receveur local i
! nrecloc                                : nombre de receveurs locaux
! nspec_acoustic                         : nombre local d'elements spectraux acoustiques

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
                                nspec_bottom, &
                                nspec_left, &
                                nspec_right, &
                                nspec_top, &
                                numabs, abs_boundary_ij, &
                                abs_boundary_normal, &
                                abs_boundary_jacobian1Dw, &
                                nelemabs, &
                                cote_abs, &
                                ib_bottom, &
                                ib_left, &
                                ib_right, &
                                ib_top, &
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
                                gammar_store_loc)


!!! Parametres fournis

! rmass_inverse_acoustic                 : matrice acoustique inversee de taille nglob
!                                          (nglob_acoustic = nglob s'il existe des elements acoustiques)
! num_phase_ispec_acoustic               : max entre nb d'element spectraux acoustiques interieur et exterieur
! phase_ispec_inner_acoustic(i,j)        : i eme element spectral acoustique interieur si j=2 exterieur si j=1
! acoustic(i)                            : vrai si l'element spectral i est acoustique
! nelem_acoustic_surface                 : nombre d'elements spectraux situes sur une surface libre acoustique
! free_ac_ispec                          : numero d'element spectral du i eme element acoustique sur surface libre
! free_surface_ij(i,j,ispec)             : i eme coordonnee du j eme point GLL de l'element spectral libre ispec
! b_reclen_potential                     : place en octet prise par b_nelem_acoustic_surface * GLLX
! any_elastic                            : vrai s'il existe des elements elastiques
! num_fluid_solid_edges                  : nombre d'elements spectraux sur une frontiere elastique/acoustique
! coupling_ac_el_ispec                   : tableau des elements spectraux frontiere ACOUSTIQUE
! coupling_ac_el_ij                      : coordonnees locales des points GLL sur la frontiere elastique/acoustique
! coupling_ac_el_normal(i,j,ispec)       : i eme coordonne de la normale au point GLL j de l'element frontiere ispec
! coupling_ac_el_jacobian1Dw(i,ispec)    : jacobienne ponderee du i eme point GLL de l'element frontiere ispec


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

!!! Parametres fournis

! rmass_inverse_elastic          : matrice elastique inversee de taille nglob_acoustic
!                                 (nglob_acoustic = nglob s'il existe des elements acoustiques)
! num_phase_ispec_elastic        : max entre nb d'element spectraux elastiques interieur et exterieur
! phase_ispec_inner_elastic(i,j) : i eme element spectral elastique interieur si j=2 exterieur si j=1
! elastic(i)                     : vrai si l'element spectral i est elastique

  ! prepares fields on GPU for elastic simulations
  !?!? JC JC here we will need to add GPU support for the new C-PML routines
  if (any_elastic) then
    ! temporary mass matrices
    allocate(rmassx(nglob_elastic),rmassz(nglob_elastic))
    rmassx(:) = rmass_inverse_elastic(1,:)
    rmassz(:) = rmass_inverse_elastic(2,:)

    call prepare_fields_elastic_device(Mesh_pointer, &
                                       rmassx,rmassz, &
                                       rho_vp,rho_vs, &
                                       num_phase_ispec_elastic,phase_ispec_inner_elastic, &
                                       ispec_is_elastic, &
                                       nspec_left, &
                                       nspec_right, &
                                       nspec_top, &
                                       nspec_bottom, &
                                       ANY_ANISOTROPY, &
                                       c11store,c12store,c13store, &
                                       c15store,c23store, &
                                       c25store,c33store,c35store,c55store, &
                                       ninterface_elastic,inum_interfaces_elastic,ATTENUATION_VISCOELASTIC, &
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

  ! synchronizes processes
  call synchronize_all()

  ! outputs GPU usage to files for all processes
  call output_free_device_memory(myrank)

  ! outputs usage for main process
  if (myrank == 0) then
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

  use constants, only: IMAIN,IEDGE1,IEDGE2,IEDGE3,IEDGE4,IBOTTOM,IRIGHT,ITOP,ILEFT

  implicit none

  ! local parameters
  integer :: i_spec_free, ipoint1D, i, j, k, ispec, ispecabs, i_source
  integer :: ispec_acoustic,ispec_elastic,iedge_acoustic,iedge_elastic
  integer :: ier,inum
  real(kind=CUSTOM_REAL) :: zxi,xgamma,jacobian1D
  real(kind=CUSTOM_REAL) :: xxi,zgamma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Initialisation variables pour routine prepare_constants_device
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

  ! checks
  if (nelemabs < 0) then
    if (myrank == 0) write(IMAIN,*) '  Warning: reading in negative nelemabs ',nelemabs,'...resetting to zero!'
    nelemabs = 0
  endif

  allocate(abs_boundary_ij(2,NGLLX,nelemabs), &
           abs_boundary_jacobian1Dw(NGLLX,nelemabs), &
           abs_boundary_normal(NDIM,NGLLX,nelemabs), &
           cote_abs(nelemabs),stat=ier)
  if (ier /= 0 ) call stop_the_code('error allocating array abs_boundary_ispec etc.')

  if (STACEY_ABSORBING_CONDITIONS) then

    do ispecabs = 1,nelemabs
      ispec = numabs(ispecabs)

      !--- left absorbing boundary
      if (codeabs(IEDGE4,ispecabs)) then
        i = 1
        do j = 1,NGLLZ
          abs_boundary_ij(1,j,ispecabs) = i
          abs_boundary_ij(2,j,ispecabs) = j

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)

          abs_boundary_normal(1,j,ispecabs) = - zgamma / jacobian1D
          abs_boundary_normal(2,j,ispecabs) = + xgamma / jacobian1D

          abs_boundary_jacobian1Dw(j,ispecabs) = jacobian1D * wzgll(j)

          cote_abs(ispecabs) = 4
        enddo
        if (ibegin_edge4(ispecabs) == 2) abs_boundary_ij(2,1,ispecabs) = 6
        if (iend_edge4(ispecabs) == 4) abs_boundary_ij(2,5,ispecabs) = 6

      !--- right absorbing boundary
      else if (codeabs(IEDGE2,ispecabs)) then
        i = NGLLX
        do j = 1,NGLLZ
          abs_boundary_ij(1,j,ispecabs) = i
          abs_boundary_ij(2,j,ispecabs) = j

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)

          abs_boundary_normal(1,j,ispecabs) = + zgamma / jacobian1D
          abs_boundary_normal(2,j,ispecabs) = - xgamma / jacobian1D

          abs_boundary_jacobian1Dw(j,ispecabs) = jacobian1D * wzgll(j)

          cote_abs(ispecabs) = 2
        enddo
        if (ibegin_edge2(ispecabs) == 2) abs_boundary_ij(2,1,ispecabs) = 6
        if (iend_edge2(ispecabs) == 4) abs_boundary_ij(2,5,ispecabs) = 6

      !--- bottom absorbing boundary
      else if (codeabs(IEDGE1,ispecabs)) then
        j = 1
        do i = 1,NGLLX
          abs_boundary_ij(1,i,ispecabs) = i
          abs_boundary_ij(2,i,ispecabs) = j

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)

          abs_boundary_normal(1,i,ispecabs) = + zxi / jacobian1D
          abs_boundary_normal(2,i,ispecabs) = - xxi / jacobian1D

          abs_boundary_jacobian1Dw(i,ispecabs) = jacobian1D * wxgll(i)

          cote_abs(ispecabs) = 1

        enddo
        if (ibegin_edge1(ispecabs) == 2 .or. codeabs_corner(1,ispecabs)) abs_boundary_ij(1,1,ispecabs) = 6
        if (iend_edge1(ispecabs) == 4 .or. codeabs_corner(2,ispecabs)) abs_boundary_ij(1,5,ispecabs) = 6

      !--- top absorbing boundary
      else if (codeabs(IEDGE3,ispecabs)) then
        j = NGLLZ
        do i = 1,NGLLX
          abs_boundary_ij(1,i,ispecabs) = i
          abs_boundary_ij(2,i,ispecabs) = j

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)

          abs_boundary_normal(1,i,ispecabs) = - zxi / jacobian1D
          abs_boundary_normal(2,i,ispecabs) = + xxi / jacobian1D

          abs_boundary_jacobian1Dw(i,ispecabs) = jacobian1D * wxgll(i)

          cote_abs(ispecabs) = 3
        enddo
        if (ibegin_edge3(ispecabs) == 2 .or. codeabs_corner(3,ispecabs)) abs_boundary_ij(1,1,ispecabs) = 6
        if (iend_edge3(ispecabs) == 4 .or. codeabs_corner(4,ispecabs)) abs_boundary_ij(1,5,ispecabs) = 6

      endif
    enddo
  endif ! STACEY_ABSORBING_CONDITIONS

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
      if (j > nsources_local) call stop_the_code('Error with the number of local sources')
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialisation parametres pour simulation elastique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! anisotropy
  ANY_ANISOTROPY = .false.
  do ispec = 1, nspec
    if (ispec_is_anisotropic(ispec) ) ANY_ANISOTROPY = .true.
  enddo

  if (ANY_ANISOTROPY) then
    ! user output
    if (myrank == 0) write(IMAIN,*) '  setting up anisotropic arrays'

    allocate(c11store(NGLLX,NGLLZ,NSPEC))
    allocate(c13store(NGLLX,NGLLZ,NSPEC))
    allocate(c15store(NGLLX,NGLLZ,NSPEC))
    allocate(c33store(NGLLX,NGLLZ,NSPEC))
    allocate(c35store(NGLLX,NGLLZ,NSPEC))
    allocate(c55store(NGLLX,NGLLZ,NSPEC))
    allocate(c12store(NGLLX,NGLLZ,NSPEC))
    allocate(c23store(NGLLX,NGLLZ,NSPEC))
    allocate(c25store(NGLLX,NGLLZ,NSPEC))

    if (assign_external_model) then
      do ispec = 1,nspec
        do j = 1,NGLLZ
          do i = 1,NGLLX
            c11store(i,j,ispec) = c11ext(i,j,ispec)
            c13store(i,j,ispec) = c13ext(i,j,ispec)
            c15store(i,j,ispec) = c15ext(i,j,ispec)
            c33store(i,j,ispec) = c33ext(i,j,ispec)
            c35store(i,j,ispec) = c35ext(i,j,ispec)
            c55store(i,j,ispec) = c55ext(i,j,ispec)
            c12store(i,j,ispec) = c12ext(i,j,ispec)
            c23store(i,j,ispec) = c23ext(i,j,ispec)
            c25store(i,j,ispec) = c25ext(i,j,ispec)
          enddo
       enddo
      enddo
    else
      do ispec = 1,nspec
        do j = 1,NGLLZ
          do i = 1,NGLLX
            c11store(i,j,ispec) = sngl(anisotropy(1,kmato(ispec)))
            c13store(i,j,ispec) = sngl(anisotropy(2,kmato(ispec)))
            c15store(i,j,ispec) = sngl(anisotropy(3,kmato(ispec)))
            c33store(i,j,ispec) = sngl(anisotropy(4,kmato(ispec)))
            c35store(i,j,ispec) = sngl(anisotropy(5,kmato(ispec)))
            c55store(i,j,ispec) = sngl(anisotropy(6,kmato(ispec)))
            c12store(i,j,ispec) = sngl(anisotropy(7,kmato(ispec)))
            c23store(i,j,ispec) = sngl(anisotropy(8,kmato(ispec)))
            c25store(i,j,ispec) = sngl(anisotropy(9,kmato(ispec)))
          enddo
        enddo
      enddo
    endif
  else
    ! dummy allocations
    allocate(c11store(1,1,1))
    allocate(c13store(1,1,1))
    allocate(c15store(1,1,1))
    allocate(c33store(1,1,1))
    allocate(c35store(1,1,1))
    allocate(c55store(1,1,1))
    allocate(c12store(1,1,1))
    allocate(c23store(1,1,1))
    allocate(c25store(1,1,1))
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  initialization done successfully'
    call flush_IMAIN()
  endif

  end subroutine init_host_to_dev_variable

!----------------------------------------------------------------------

  end subroutine prepare_GPU
