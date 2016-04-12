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
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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


  subroutine prepare_timerun_GPU()

  use specfem_par
  use specfem_par_gpu

  implicit none

  ! local parameters
  real :: free_mb,used_mb,total_mb
  integer :: nspec_elastic
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: cosrot_irecf, sinrot_irecf

  ! GPU_MODE now defined in Par_file
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) "GPU Preparing Fields and Constants on Device."
    call flush_IMAIN()
  endif

  ! initializes arrays
  call init_host_to_dev_variable()

  ! number of purely elastic elements
  nspec_elastic = nspec - count_nspec_acoustic

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
! recloc(i)                              : convertisseur numero rec local i => numero rec global
! ispec_selected_rec(i)                  : numero d'element spectral du receveur i
! nrecloc                                : nombre de receveurs locaux
! count_nspec_acoustic                   : nombre local d'elements spectraux acoustiques

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
                                wxgll,&
                                STACEY_BOUNDARY_CONDITIONS, &
                                nspec_bottom,&
                                nspec_left,&
                                nspec_right,&
                                nspec_top,&
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
                                sourcearray_loc,source_time_function_loc,&
                                NSTEP,ispec_selected_source_loc, &
                                recloc, ispec_selected_rec, &
                                nrec, nrecloc, &
                                cosrot_irecf,sinrot_irecf,&
                                SIMULATION_TYPE, &
                                USE_MESH_COLORING_GPU, &
                                count_nspec_acoustic,nspec_elastic,&
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
! b_reclen                               : place en octet prise par b_nelem_acoustic_surface * GLLX
! any_elastic                            : vrai s'il existe des elements elastiques
! num_fluid_solid_edges                  : nombre d'elements spectraux sur une frontiere elastique/acoustique
! coupling_ac_el_ispec                   : tableau des elements spectraux frontiere ACOUSTIQUE
! coupling_ac_el_ij                      : coordonnees locales des points GLL sur la frontiere elastique/acoustique
! coupling_ac_el_normal(i,j,ispec)       : i eme coordonne de la normale au point GLL j de l'element frontiere ispec
! coupling_ac_el_jacobian1Dw(i,ispec)    : jacobienne ponderee du i eme point GLL de l'element frontiere ispec
! num_colors_outer_acoustic              : a initialiser plus tard quand USE_COLOR_MESH sera implemente
! num_colors_inner_acoustic              : a initialiser plus tard quand USE_COLOR_MESH sera implemente
! num_elem_colors_acoustic               : a initialiser plus tard quand USE_COLOR_MESH sera implemente


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
                                ninterface_acoustic,inum_interfaces_acoustic, &
                                num_colors_outer_acoustic,num_colors_inner_acoustic, &
                                num_elem_colors_acoustic)

    if (SIMULATION_TYPE == 3) then
      ! safety check
      if (APPROXIMATE_HESS_KL) then
        stop 'Sorry, approximate acoustic hessian kernels not yet fully implemented for GPU simulations!'
      endif
      call prepare_fields_acoustic_adj_dev(Mesh_pointer,APPROXIMATE_HESS_KL)
    endif
  endif

!!! Parametres fournis

! rmass_inverse_elastic          : matrice elastique inversee de taille nglob_acoustic
!                                 (nglob_acoustic = nglob s'il existe des elements acoustiques)
! num_phase_ispec_elastic        : max entre nb d'element spectraux elastiques interieur et exterieur
! phase_ispec_inner_elastic(i,j) : i eme element spectral elastique interieur si j=2 exterieur si j=1
! elastic(i)                     : vrai si l'element spectral i est elastique
! num_colors_outer_elastic       : a initialiser plus tard quand USE_COLOR_MESH sera implemente
! num_colors_inner_elastic       : a initialiser plus tard quand USE_COLOR_MESH sera implemente
! num_elem_colors_elastic        : a initialiser plus tard quand USE_COLOR_MESH sera implemente

  ! prepares fields on GPU for elastic simulations
  !?!? JC JC here we will need to add GPU support for the new C-PML routines
  if (any_elastic) then
    call prepare_fields_elastic_device(Mesh_pointer, &
                                rmass_inverse_elastic_one,rmass_inverse_elastic_three, &
                                rho_vp,rho_vs, &
                                num_phase_ispec_elastic,phase_ispec_inner_elastic, &
                                ispec_is_elastic, &
                                nspec_left,&
                                nspec_right,&
                                nspec_top,&
                                nspec_bottom,&
                                any_acoustic, &
                                num_colors_outer_elastic,num_colors_inner_elastic, &
                                num_elem_colors_elastic, &
                                ANY_ANISOTROPY, &
                                c11store,c12store,c13store, &
                                c15store,c23store, &
                                c25store,c33store,c35store,c55store,ninterface_elastic,inum_interfaces_elastic)


    if (SIMULATION_TYPE == 3) then
      ! safety check
      if (APPROXIMATE_HESS_KL) then
        stop 'Sorry, approximate elastic hessian kernels not yet fully implemented for GPU simulations!'
      endif
      call prepare_fields_elastic_adj_dev(Mesh_pointer,NDIM*NGLOB_AB,APPROXIMATE_HESS_KL)
    endif
  endif

  ! prepares fields on GPU for poroelastic simulations
  if (any_poroelastic) then
    stop 'todo poroelastic simulations on GPU'
  endif


  ! prepares needed receiver array for adjoint runs
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) &
    call prepare_sim2_or_3_const_device(Mesh_pointer, &
                                which_proc_receiver,nrecloc,nrec,source_adjointe,NSTEP)


  ! synchronizes processes
  call synchronize_all()


  ! puts acoustic initial fields onto GPU
  if (any_acoustic) then
    call transfer_fields_ac_to_device(NGLOB_AB,minus_int_int_pressure_acoustic, &
                                      minus_int_pressure_acoustic,minus_pressure_acoustic,Mesh_pointer)



    if (SIMULATION_TYPE == 3) &
      call transfer_b_fields_ac_to_device(NGLOB_AB,b_minus_int_int_pressure_acoustic, &
                                          b_minus_int_pressure_acoustic,b_minus_pressure_acoustic,Mesh_pointer)
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

  ! allocates arrays for mpi transfers
  allocate(tab_requests_send_recv_scalar(2*ninterface))
  allocate(b_tab_requests_send_recv_scalar(2*ninterface))
  allocate(tab_requests_send_recv_vector(2*ninterface))
  allocate(b_tab_requests_send_recv_vector(2*ninterface))

  allocate(buffer_send_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,ninterface))
  allocate(b_buffer_send_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,ninterface))
  allocate(buffer_recv_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,ninterface))
  allocate(b_buffer_recv_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,ninterface))
  allocate(buffer_send_vector_ext_mesh(2,max_nibool_interfaces_ext_mesh,ninterface))
  allocate(b_buffer_send_vector_ext_mesh(2,max_nibool_interfaces_ext_mesh,ninterface))
  allocate(buffer_recv_vector_ext_mesh(2,max_nibool_interfaces_ext_mesh,ninterface))
  allocate(b_buffer_recv_vector_ext_mesh(2,max_nibool_interfaces_ext_mesh,ninterface))

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

  implicit none

  ! local parameters
  integer :: i_spec_free, ipoint1D, i, j, k, ispec, ispecabs, i_source, ispec_inner, ispec_outer
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

  ! converts to custom_real
  deltatf = sngl(deltat)
  deltatover2f = sngl(deltatover2)
  deltatsquareover2f = sngl(deltatsquareover2)

  b_deltatf = sngl(b_deltat)
  b_deltatover2f = sngl(b_deltatover2)
  b_deltatsquareover2f = sngl(b_deltatsquareover2)

  NSPEC_AB = nspec
  NGLOB_AB = nglob

  ! user output
  if (myrank == 0) write(IMAIN,*) '  number of mpi interfaces = ',ninterface

  ! mpi interfaces
  allocate(ibool_interfaces_ext_mesh(max_nibool_interfaces_ext_mesh,ninterface))
  ibool_interfaces_ext_mesh(:,:) = 0
  do j = 1,ninterface
    do i = 1,nibool_interfaces_ext_mesh(j)
      ibool_interfaces_ext_mesh(i,j)=ibool_interfaces_ext_mesh_init(i,j)
    enddo
  enddo

  ! user output
  if (myrank == 0) write(IMAIN,*) '  number of acoustic elements at free surface = ',nelem_acoustic_surface

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
           abs_boundary_normal(NDIM,NGLLX,nelemabs),&
           cote_abs(nelemabs),stat=ier)
  if (ier /= 0 ) stop 'error allocating array abs_boundary_ispec etc.'

  if(STACEY_BOUNDARY_CONDITIONS) then

    do ispecabs = 1,nelemabs
      ispec = numabs(ispecabs)

      !--- left absorbing boundary
      if(codeabs(IEDGE4,ispecabs)) then
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
        if (ibegin_edge4(ispecabs)==2) abs_boundary_ij(2,1,ispecabs) = 6
        if (iend_edge4(ispecabs)==4) abs_boundary_ij(2,5,ispecabs) = 6

      !--- right absorbing boundary
      else if(codeabs(IEDGE2,ispecabs)) then
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
        if (ibegin_edge2(ispecabs)==2) abs_boundary_ij(2,1,ispecabs) = 6
        if (iend_edge2(ispecabs)==4) abs_boundary_ij(2,5,ispecabs) = 6

      !--- bottom absorbing boundary
      else if(codeabs(IEDGE1,ispecabs)) then
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
        if (ibegin_edge1(ispecabs)==2  .or. codeabs_corner(1,ispecabs)) abs_boundary_ij(1,1,ispecabs) = 6
        if (iend_edge1(ispecabs)==4 .or. codeabs_corner(2,ispecabs))    abs_boundary_ij(1,5,ispecabs) = 6

      !--- top absorbing boundary
      else if(codeabs(IEDGE3,ispecabs)) then
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
        if (ibegin_edge3(ispecabs)==2 .or. codeabs_corner(3,ispecabs)) abs_boundary_ij(1,1,ispecabs) = 6
        if (iend_edge3(ispecabs)==4 .or. codeabs_corner(4,ispecabs))   abs_boundary_ij(1,5,ispecabs) = 6

      endif
    enddo
  endif ! STACEY_BOUNDARY_CONDITIONS

  ! user output
  if (myrank == 0) write(IMAIN,*) '  number of sources = ',NSOURCES

  ! counts sources in this process slice
  nsources_local = 0
  do i = 1, NSOURCES
    if (is_proc_source(i) == 1) then
      nsources_local = nsources_local + 1
    endif
  enddo

  ! user output
  if (myrank == 0) write(IMAIN,*) '  number of sources in this process slice = ',nsources_local

  allocate(source_time_function_loc(nsources_local,NSTEP))
  allocate(ispec_selected_source_loc(nsources_local))
  j = 0
  do i = 1, NSOURCES
    if (is_proc_source(i) == 1) then
      if (j>nsources_local) stop 'Error with the number of local sources'
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
    if (is_proc_source(i_source) == 1) then
      ! source belongs to this process
      k = k + 1
      if (source_type(i_source) == 1) then
        if (ispec_is_acoustic(ispec_selected_source(i_source))) then
          ! acoustic domain
          do j = 1,NGLLZ
            do i = 1,NGLLX
              sourcearray_loc(k,1,i,j) = sngl(hxis_store(i_source,i) * hgammas_store(i_source,j))
            enddo
          enddo
        else if (ispec_is_elastic(ispec_selected_source(i_source))) then
          ! elastic domain
          do j = 1,NGLLZ
            do i = 1,NGLLX
              sourcearray_loc(k,1,i,j) = - sngl(sin(anglesource(i_source)) * hxis_store(i_source,i) * hgammas_store(i_source,j))
              sourcearray_loc(k,2,i,j) = sngl(cos(anglesource(i_source)) * hxis_store(i_source,i) * hgammas_store(i_source,j))
            enddo
          enddo
        endif
      else
        sourcearray_loc(k,:,:,:) = sourcearray(i_source,:,:,:)
        sourcearray_loc(k,:,:,:) = sourcearray(i_source,:,:,:)
      endif ! Source_type
    endif ! is_proc_source
  enddo

  ! converts to custom_real arrays
  allocate(cosrot_irecf(nrecloc), &
           sinrot_irecf(nrecloc))
  do i = 1,nrecloc
    cosrot_irecf(i) = sngl(cosrot_irec(i))
    sinrot_irecf(i) = sngl(sinrot_irec(i))
  enddo

  ! determines inner elements
  allocate(ispec_is_inner(nspec))
  ispec_is_inner(:) = .false.

  ! loop over spectral elements
  do ispec_inner = 1,nspec_inner
    ! get global numbering for inner or outer elements
    ispec = ispec_inner_to_glob(ispec_inner)
    ispec_is_inner(ispec) = .true.
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Init pour prepare acoustique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! sets up elements for loops in acoustic simulations
  nspec_inner_acoustic = 0
  nspec_outer_acoustic = 0
  if (any_acoustic) then
    ! user output
    if (myrank == 0) write(IMAIN,*) '  determining phases for acoustic domain'

    ! counts inner and outer elements
    do ispec = 1, nspec
      if (ispec_is_acoustic(ispec)) then
        if (ispec_is_inner(ispec) .eqv. .true.) then
          nspec_inner_acoustic = nspec_inner_acoustic + 1
        else
          nspec_outer_acoustic = nspec_outer_acoustic + 1
        endif
      endif
    enddo

    ! stores indices of inner and outer elements for faster(?) computation
    num_phase_ispec_acoustic = max(nspec_inner_acoustic,nspec_outer_acoustic)
    if (num_phase_ispec_acoustic < 0 ) stop 'Error acoustic simulation: num_phase_ispec_acoustic is < zero'

    allocate( phase_ispec_inner_acoustic(num_phase_ispec_acoustic,2),stat=ier)
    if (ier /= 0 ) stop 'Error allocating array phase_ispec_inner_acoustic'
    phase_ispec_inner_acoustic(:,:) = 0

    ispec_inner = 0
    ispec_outer = 0
    do ispec = 1, nspec
      if (ispec_is_acoustic(ispec)) then
        if (ispec_is_inner(ispec) .eqv. .true.) then
          ispec_inner = ispec_inner + 1
          phase_ispec_inner_acoustic(ispec_inner,2) = ispec
        else
          ispec_outer = ispec_outer + 1
          phase_ispec_inner_acoustic(ispec_outer,1) = ispec
        endif
      endif
    enddo
  else
    ! allocates dummy array
    num_phase_ispec_acoustic = 0
    allocate( phase_ispec_inner_acoustic(num_phase_ispec_acoustic,2),stat=ier)
    if (ier /= 0 ) stop 'Error allocating dummy array phase_ispec_inner_acoustic'
    phase_ispec_inner_acoustic(:,:) = 0
  endif

  ! free surface
  allocate(free_surface_ij(2,NGLLX,nelem_acoustic_surface))

  do i_spec_free = 1, nelem_acoustic_surface
    if (acoustic_surface(2,i_spec_free) == acoustic_surface(3,i_spec_free)) then
      do j =1,5
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
      do j =1,5
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

      else if (iedge_acoustic ==ILEFT) then
        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        coupling_ac_el_normal(1,ipoint1D,inum) = - zgamma / jacobian1D
        coupling_ac_el_normal(2,ipoint1D,inum) = + xgamma / jacobian1D
        coupling_ac_el_jacobian1Dw(ipoint1D,inum) = jacobian1D * wzgll(j)

      else if (iedge_acoustic ==IRIGHT) then
        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        coupling_ac_el_normal(1,ipoint1D,inum) = + zgamma / jacobian1D
        coupling_ac_el_normal(2,ipoint1D,inum) = - xgamma / jacobian1D
        coupling_ac_el_jacobian1Dw(ipoint1D,inum) = jacobian1D * wzgll(j)
      endif
    enddo
  enddo

  ! coloring (dummy)
  num_colors_outer_acoustic = 0
  num_colors_inner_acoustic = 0
  allocate(num_elem_colors_acoustic(1))
  num_elem_colors_acoustic(1) = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialisation parametres pour simulation elastique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! sets up elements for loops in elastic simulations
  nspec_inner_elastic = 0
  nspec_outer_elastic = 0
  if (any_elastic) then
    ! user output
    if (myrank == 0) write(IMAIN,*) '  determining phases for elastic domain'

    ! counts inner and outer elements
    do ispec = 1, nspec
      if (ispec_is_elastic(ispec)) then
        if (ispec_is_inner(ispec) .eqv. .true.) then
          nspec_inner_elastic = nspec_inner_elastic + 1
        else
          nspec_outer_elastic = nspec_outer_elastic + 1
        endif
      endif
    enddo

    ! stores indices of inner and outer elements for faster(?) computation
    num_phase_ispec_elastic = max(nspec_inner_elastic,nspec_outer_acoustic)
    if (num_phase_ispec_elastic < 0 ) stop 'error elastic simulation: num_phase_ispec_elastic is < zero'

    allocate( phase_ispec_inner_elastic(num_phase_ispec_elastic,2),stat=ier)
    if (ier /= 0 ) stop 'Error allocating array phase_ispec_inner_elastic'
    phase_ispec_inner_elastic(:,:) = 0

    ispec_inner = 0
    ispec_outer = 0
    do ispec = 1, nspec
      if (ispec_is_elastic(ispec)) then
        if (ispec_is_inner(ispec) .eqv. .true.) then
          ispec_inner = ispec_inner + 1
          phase_ispec_inner_elastic(ispec_inner,2) = ispec
        else
          ispec_outer = ispec_outer + 1
          phase_ispec_inner_elastic(ispec_outer,1) = ispec
        endif
      endif
    enddo
  else
    ! allocates dummy array
    num_phase_ispec_elastic = 0
    allocate( phase_ispec_inner_elastic(num_phase_ispec_elastic,2),stat=ier)
    if (ier /= 0 ) stop 'Error allocating dummy array phase_ispec_inner_elastic'
    phase_ispec_inner_elastic(:,:) = 0
  endif

  ! coloring (dummy)
  num_colors_outer_elastic = 0
  num_colors_inner_elastic = 0
  allocate(num_elem_colors_elastic(1))
  num_elem_colors_elastic(1)=0

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

  end subroutine prepare_timerun_GPU
