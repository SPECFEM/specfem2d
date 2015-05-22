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

 subroutine prepare_timerun_GPU()

  use specfem_par

  implicit none

  ! local parameters
  real :: free_mb,used_mb,total_mb
  integer :: nspec_elastic

  nspec_elastic = nspec - count_nspec_acoustic

  ! GPU_MODE now defined in Par_file
  if(myrank == 0 ) then
    write(IOUT,*)
    write(IOUT,*) "GPU Preparing Fields and Constants on Device."
    call flush_IOUT()
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

! rmass_inverse_acoustic                 : matrice acoustique inversee de taille nglob (nglob_acoustic = nglob s'il existe des elements acoustiques)
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
! num_colors_outer_acoustic              : a initialiser plus tard quand USE_COLOR_MESH sera implemente
! num_colors_inner_acoustic              : a initialiser plus tard quand USE_COLOR_MESH sera implemente
! num_elem_colors_acoustic               : a initialiser plus tard quand USE_COLOR_MESH sera implemente


  ! prepares fields on GPU for acoustic simulations
  if( any_acoustic ) then
    call prepare_fields_acoustic_device(Mesh_pointer, &
                                rmass_inverse_acoustic,rhostore,kappastore, &
                                num_phase_ispec_acoustic,phase_ispec_inner_acoustic, &
                                acoustic, &
                                nelem_acoustic_surface, &
                                free_ac_ispec,free_surface_ij, &
                                any_elastic, num_fluid_solid_edges, &
                                coupling_ac_el_ispec,coupling_ac_el_ij, &
                                coupling_ac_el_normal,coupling_ac_el_jacobian1Dw, &
                                ninterface_acoustic,inum_interfaces_acoustic, &
                                num_colors_outer_acoustic,num_colors_inner_acoustic, &
                                num_elem_colors_acoustic)

    if( SIMULATION_TYPE == 3 ) &
      call prepare_fields_acoustic_adj_dev(Mesh_pointer, &
                                APPROXIMATE_HESS_KL)

  endif

!!! Parametres fournis

! rmass_inverse_elastic                 : matrice elastique inversee de taille nglob_acoustic (nglob_acoustic = nglob s'il existe des elements acoustiques)
! num_phase_ispec_elastic               : max entre nb d'element spectraux elastiques interieur et exterieur
! phase_ispec_inner_elastic(i,j)        : i eme element spectral elastique interieur si j=2 exterieur si j=1
! elastic(i)                            : vrai si l'element spectral i est elastique
! num_colors_outer_elastic              : a initialiser plus tard quand USE_COLOR_MESH sera implemente
! num_colors_inner_elastic              : a initialiser plus tard quand USE_COLOR_MESH sera implemente
! num_elem_colors_elastic               : a initialiser plus tard quand USE_COLOR_MESH sera implemente

  ! prepares fields on GPU for elastic simulations
  !?!? JC JC here we will need to add GPU support for the new C-PML routines
  if( any_elastic ) then
    call prepare_fields_elastic_device(Mesh_pointer, &
                                rmass_inverse_elastic_one,rmass_inverse_elastic_three, &
                                rho_vp,rho_vs, &
                                num_phase_ispec_elastic,phase_ispec_inner_elastic, &
                                elastic, &
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


    if( SIMULATION_TYPE == 3 ) &
      call prepare_fields_elastic_adj_dev(Mesh_pointer, &
                                NDIM*NGLOB_AB, &
                                APPROXIMATE_HESS_KL)

  endif

  ! prepares fields on GPU for poroelastic simulations
  if( any_poroelastic ) then
    stop 'todo poroelastic simulations on GPU'
  endif


  ! prepares needed receiver array for adjoint runs
  if( SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3 ) &
    call prepare_sim2_or_3_const_device(Mesh_pointer, &
                                which_proc_receiver,nrecloc,nrec,source_adjointe,NSTEP)


  ! synchronizes processes
  call sync_all()


  ! puts acoustic initial fields onto GPU
  if( any_acoustic ) then
    call transfer_fields_ac_to_device(NGLOB_AB,potential_acoustic, &
                                      potential_dot_acoustic,potential_dot_dot_acoustic,Mesh_pointer)



    if( SIMULATION_TYPE == 3 ) &
      call transfer_b_fields_ac_to_device(NGLOB_AB,b_potential_acoustic, &
                                          b_potential_dot_acoustic,b_potential_dot_dot_acoustic,Mesh_pointer)
  endif

  ! puts elastic initial fields onto GPU
  if( any_elastic ) then
    ! transfer forward and backward fields to device with initial values
    call transfer_fields_el_to_device(NDIM*NGLOB_AB,displ_2D,veloc_2D,accel_2D,Mesh_pointer)
    if(SIMULATION_TYPE == 3) &
      call transfer_b_fields_to_device(NDIM*NGLOB_AB,b_displ_2D,b_veloc_2D,b_accel_2D,Mesh_pointer)
  endif

  ! synchronizes processes
  call sync_all()

  ! outputs GPU usage to files for all processes
  call output_free_device_memory(myrank)

  ! outputs usage for main process
  if( myrank == 0 ) then
    call get_free_device_memory(free_mb,used_mb,total_mb)
    write(IOUT,*) "GPU usage: free  =",free_mb," MB",nint(free_mb/total_mb*100.0),"%"
    write(IOUT,*) "           used  =",used_mb," MB",nint(used_mb/total_mb*100.0),"%"
    write(IOUT,*) "           total =",total_mb," MB",nint(total_mb/total_mb*100.0),"%"
    write(IOUT,*)
    call flush_IOUT()
  endif

  end subroutine prepare_timerun_GPU
