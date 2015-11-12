
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

subroutine compute_forces_viscoelastic(accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old, &
                                       x0_source, z0_source,f0,v0x_left,v0z_left,v0x_right,v0z_right,&
                                       v0x_bot,v0z_bot,t0x_left,t0z_left,t0x_right,t0z_right,t0x_bot,t0z_bot,&
                                       nleft,nright,nbot,PML_BOUNDARY_CONDITIONS,e1,e11,e13)

  ! compute forces for the elastic elements

  use specfem_par, only: p_sv,nglob,nspec,nelemabs,it,anyabs,assign_external_model, &
                         initialfield,ATTENUATION_VISCOELASTIC_SOLID,nspec_allocate,N_SLS,anglesource, &
                         ibool,kmato,numabs,elastic,codeabs,codeabs_corner, &
                         density,poroelastcoef,xix,xiz,gammax,gammaz, &
                         jacobian,vpext,vsext,rhoext,c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,c22ext, &
                         anisotropic,anisotropy, &
                         e1_LDDRK,e11_LDDRK,e13_LDDRK,alpha_LDDRK,beta_LDDRK,c_LDDRK, &
                         e1_initial_rk,e11_initial_rk,e13_initial_rk,e1_force_RK, e11_force_RK, e13_force_RK, &
                         hprime_xx,hprimewgll_xx,hprime_zz,hprimewgll_zz,wxgll,wzgll, &
                         AXISYM,is_on_the_axis,hprimeBar_xx,hprimeBarwglj_xx,xiglj,wxglj, &
                         inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2,N_SLS, &
                         deltat,coord,add_Bielak_conditions, &
                         A_plane, B_plane, C_plane, anglesource_refl, c_inc, c_refl, time_offset, &
                         over_critical_angle,SIMULATION_TYPE,SAVE_FORWARD,b_absorb_elastic_left,&
                         b_absorb_elastic_right,b_absorb_elastic_bottom,b_absorb_elastic_top,&
                         ib_left,ib_right,ib_bottom,ib_top,&
                         stage_time_scheme,i_stage,ADD_SPRING_TO_STACEY,x_center_spring,z_center_spring,&
                         is_PML,nspec_PML,spec_to_PML,region_CPML,rmemory_duz_dz_LDDRK, &
                         K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store, &
                         rmemory_displ_elastic,rmemory_dux_dx,rmemory_dux_dz,rmemory_duz_dx,rmemory_duz_dz, &
                         rmemory_dux_dx_prime,rmemory_dux_dz_prime,rmemory_duz_dx_prime,rmemory_duz_dz_prime, &
                         rmemory_displ_elastic_LDDRK,rmemory_dux_dx_LDDRK,rmemory_dux_dz_LDDRK,rmemory_duz_dx_LDDRK,&
                         ROTATE_PML_ACTIVATE,ROTATE_PML_ANGLE,STACEY_BOUNDARY_CONDITIONS,acoustic,time_stepping_scheme

  implicit none
  include "constants.h"

  real(kind=CUSTOM_REAL), dimension(3,nglob) :: accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_allocate,N_SLS) :: e1,e11,e13

  ! for analytical initial plane wave for Bielak's conditions
  double precision x0_source, z0_source,f0

  integer :: nleft, nright, nbot
  double precision, dimension(nleft) :: v0x_left,v0z_left,t0x_left,t0z_left
  double precision, dimension(nright) :: v0x_right,v0z_right,t0x_right,t0z_right
  double precision, dimension(nbot) :: v0x_bot,v0z_bot,t0x_bot,t0z_bot

  ! CPML coefficients and memory variables
  logical :: PML_BOUNDARY_CONDITIONS

  !---
  !--- local variables
  !---

  integer :: ispec,i,j,k,iglob,ispecabs,ibegin,iend

  ! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duy_dxi,duy_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dxl,duy_dxl,duz_dxl,dux_dzl,duy_dzl,duz_dzl
  real(kind=CUSTOM_REAL) :: dux_dxl_prime,duz_dxl_prime,dux_dzl_prime,duz_dzl_prime
  real(kind=CUSTOM_REAL) :: theta,ct,st
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_xy,sigma_xz,sigma_zy,sigma_zz,sigma_zx
  real(kind=CUSTOM_REAL) :: sigma_xx_prime,sigma_xz_prime,sigma_zz_prime,sigma_zx_prime
  real(kind=CUSTOM_REAL) :: nx,nz,vx,vy,vz,vn,rho_vp,rho_vs,tx,ty,tz,weight,xxi,zxi,xgamma,zgamma,jacobian1D
  real(kind=CUSTOM_REAL) :: displx,disply,displz,displn,spring_position,displtx,displty,displtz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2,tempy1,tempy2,tempz1,tempz2
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: tempx3
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: sigma_thetatheta
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) :: e1_sum,e11_sum,e13_sum
  integer :: i_sls

  ! nsub1 denotes discrete time step n-1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: dux_dxl_n,duz_dzl_n,duz_dxl_n,dux_dzl_n, &
                                                        dux_dxl_nsub1,duz_dzl_nsub1,duz_dxl_nsub1,dux_dzl_nsub1

  double precision :: coef0,coef1,coef2

  ! material properties of the elastic medium
  real(kind=CUSTOM_REAL) :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic, &
    lambdaplus2mu_unrelaxed_elastic,kappal,cpl,csl,rhol,lambdalplusmul_unrelaxed_elastic

  ! for attenuation
  real(kind=CUSTOM_REAL) :: phinu1,phinu2,theta_n_u,theta_nsub1_u
  double precision :: tauinvnu1,tauinvnu2

  ! for anisotropy
  double precision ::  c11,c15,c13,c33,c35,c55,c12,c23,c25,c22

  ! for analytical initial plane wave for Bielak's conditions
  double precision :: veloc_horiz,veloc_vert,dxUx,dzUx,dxUz,dzUz,traction_x_t0,traction_z_t0
  integer count_left,count_right,count_bottom

  integer :: ifirstelem,ilastelem

  ! CPML coefficients and memory variables
  integer :: ispec_PML
  real(kind=CUSTOM_REAL), dimension(3,NGLLX,NGLLZ) :: accel_elastic_PML
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) ::PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl,&
                           PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old

  real(kind=CUSTOM_REAL) :: dux_dxi_old,dux_dgamma_old,duz_dxi_old,duz_dgamma_old

  double precision :: time_n,time_nsub1
  double precision :: kappa_x,kappa_z,d_x,d_z,alpha_x,alpha_z,beta_x,beta_z, &
                      A5,A6,A7, bb_zx_1,bb_zx_2,coef0_zx_1,coef1_zx_1,coef2_zx_1,coef0_zx_2,coef1_zx_2,coef2_zx_2,&
                      A8,A9,A10,bb_xz_1,bb_xz_2,coef0_xz_1,coef1_xz_1,coef2_xz_1,coef0_xz_2,coef1_xz_2,coef2_xz_2,&
                      A0,A1,A2,A3,A4,bb_1,coef0_1,coef1_1,coef2_1,bb_2,coef0_2,coef1_2,coef2_2
  integer :: CPML_region_local,singularity_type_zx,singularity_type_xz,singularity_type

  ! temporary RK4 variable
  real(kind=CUSTOM_REAL) :: weight_rk

  !!!update memory variable in viscoelastic simulation
  if( ATTENUATION_VISCOELASTIC_SOLID ) then

    ! compute Grad(displ_elastic) at time step n for attenuation
    call compute_gradient_attenuation(displ_elastic,dux_dxl_n,duz_dxl_n, &
          dux_dzl_n,duz_dzl_n,xix,xiz,gammax,gammaz,ibool,elastic,hprime_xx,hprime_zz,nspec,nglob)

    ! compute Grad(disp_elastic_old) at time step n-1 for attenuation
    call compute_gradient_attenuation(displ_elastic_old,dux_dxl_nsub1,duz_dxl_nsub1, &
          dux_dzl_nsub1,duz_dzl_nsub1,xix,xiz,gammax,gammaz,ibool,elastic,hprime_xx,hprime_zz,nspec,nglob)

    ! loop over spectral elements
    do ispec = 1,nspec

      ! attenuation is not implemented in acoustic (i.e. fluid) media for now, only in viscoelastic (i.e. solid) media
      if( acoustic(ispec) ) cycle

      if( (.not. PML_BOUNDARY_CONDITIONS) .or. (PML_BOUNDARY_CONDITIONS .and. (.not. is_PML(ispec))) ) then
        do j=1,NGLLZ
        do i=1,NGLLX

          ! convention to indicate that Q = 9999 in that element i.e. that there is no viscoelasticity in that element
          if( inv_tau_sigma_nu1(i,j,ispec,1) < 0.) cycle

          theta_n_u = dux_dxl_n(i,j,ispec) + duz_dzl_n(i,j,ispec)
          theta_nsub1_u = dux_dxl_nsub1(i,j,ispec) + duz_dzl_nsub1(i,j,ispec)

          ! loop on all the standard linear solids
          do i_sls = 1,N_SLS
            phinu1 = phi_nu1(i,j,ispec,i_sls)
            tauinvnu1 = inv_tau_sigma_nu1(i,j,ispec,i_sls)
            phinu2 = phi_nu2(i,j,ispec,i_sls)
            tauinvnu2 = inv_tau_sigma_nu2(i,j,ispec,i_sls)

            ! update e1, e11, e13 in convolution formation with modified recursive convolution scheme on basis of
            ! second-order accurate convolution term calculation from equation (21) of
            ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
            ! Anisotropic-Medium PML for Vector FETD With Modified Basis Functions,
            ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)
            if( stage_time_scheme == 1 ) then

! update the memory variables using a convolution or using a differential equation
! From Zhinan Xie and Dimitri Komatitsch:
! For cases in which a value of tau_sigma is small, then its inverse is large,
! which may result in a in stiff ordinary differential equation to solve;
! in such a case, resorting to the convolution formulation may be better (?)
              if(CONVOLUTION_MEMORY_VARIABLES) then
                call compute_coef_convolution(tauinvnu1,deltat,coef0,coef1,coef2)

                e1(i,j,ispec,i_sls) = coef0 * e1(i,j,ispec,i_sls) + &
                                      phinu1 * (coef1 * theta_n_u + coef2 * theta_nsub1_u)

                call compute_coef_convolution(tauinvnu2,deltat,coef0,coef1,coef2)

                e11(i,j,ispec,i_sls) = coef0 * e11(i,j,ispec,i_sls) + &
                                       phinu2 * (coef1 * (dux_dxl_n(i,j,ispec)-theta_n_u/TWO) + &
                                                 coef2 * (dux_dxl_nsub1(i,j,ispec)-theta_nsub1_u/TWO))

                e13(i,j,ispec,i_sls) = coef0 * e13(i,j,ispec,i_sls) + &
                                       phinu2 * (coef1 * (dux_dzl_n(i,j,ispec) + duz_dxl_n(i,j,ispec)) + &
                                                 coef2 * (dux_dzl_nsub1(i,j,ispec) + duz_dxl_nsub1(i,j,ispec)))
              else
                e1(i,j,ispec,i_sls) = e1(i,j,ispec,i_sls) + deltat * &
                     (- e1(i,j,ispec,i_sls)*tauinvnu1 + phinu1 * theta_n_u)

                e11(i,j,ispec,i_sls) = e11(i,j,ispec,i_sls) + deltat * &
                     (- e11(i,j,ispec,i_sls)*tauinvnu2 + phinu2 * (dux_dxl_n(i,j,ispec)-theta_n_u/TWO))

                e13(i,j,ispec,i_sls) = e13(i,j,ispec,i_sls) + deltat * &
                   (- e13(i,j,ispec,i_sls)*tauinvnu2 + phinu2 * (dux_dzl_n(i,j,ispec) + duz_dxl_n(i,j,ispec)))
              endif

            endif

            ! update e1, e11, e13 in ADE formation with fourth-order LDDRK scheme
            if( stage_time_scheme == 6 ) then
              e1_LDDRK(i,j,ispec,i_sls) = alpha_LDDRK(i_stage) * e1_LDDRK(i,j,ispec,i_sls) + &
                                          deltat * (theta_n_u * phinu1 - e1(i,j,ispec,i_sls) * tauinvnu1)
              e1(i,j,ispec,i_sls) = e1(i,j,ispec,i_sls) + beta_LDDRK(i_stage) * e1_LDDRK(i,j,ispec,i_sls)

              e11_LDDRK(i,j,ispec,i_sls) = alpha_LDDRK(i_stage) * e11_LDDRK(i,j,ispec,i_sls) + &
                                           deltat * ((dux_dxl_n(i,j,ispec)-theta_n_u/TWO) * phinu2) - &
                                           deltat * (e11(i,j,ispec,i_sls) * tauinvnu2)
              e11(i,j,ispec,i_sls) = e11(i,j,ispec,i_sls)+beta_LDDRK(i_stage)*e11_LDDRK(i,j,ispec,i_sls)

              e13_LDDRK(i,j,ispec,i_sls) = alpha_LDDRK(i_stage) * e13_LDDRK(i,j,ispec,i_sls) + &
                                           deltat * ((dux_dzl_n(i,j,ispec) + duz_dxl_n(i,j,ispec))*phinu2) - &
                                           deltat * (e13(i,j,ispec,i_sls) * tauinvnu2)
              e13(i,j,ispec,i_sls) = e13(i,j,ispec,i_sls)+beta_LDDRK(i_stage) * e13_LDDRK(i,j,ispec,i_sls)
            endif

            ! update e1, e11, e13 in ADE formation with classical fourth-order Runge-Kutta scheme
            if( stage_time_scheme == 4 ) then
              e1_force_RK(i,j,ispec,i_sls,i_stage) = deltat * (theta_n_u * phinu1 - e1(i,j,ispec,i_sls) * tauinvnu1)

              if( i_stage==1 .or. i_stage==2 .or. i_stage==3 ) then
                if( i_stage == 1) weight_rk = 0.5_CUSTOM_REAL
                if( i_stage == 2) weight_rk = 0.5_CUSTOM_REAL
                if( i_stage == 3) weight_rk = 1._CUSTOM_REAL

                if( i_stage==1) e1_initial_rk(i,j,ispec,i_sls) = e1(i,j,ispec,i_sls)
                e1(i,j,ispec,i_sls) = e1_initial_rk(i,j,ispec,i_sls) + weight_rk * e1_force_RK(i,j,ispec,i_sls,i_stage)
              else if( i_stage==4 ) then
                e1(i,j,ispec,i_sls) = e1_initial_rk(i,j,ispec,i_sls) + 1._CUSTOM_REAL / 6._CUSTOM_REAL * &
                                      (e1_force_RK(i,j,ispec,i_sls,1) + 2._CUSTOM_REAL * e1_force_RK(i,j,ispec,i_sls,2) + &
                                       2._CUSTOM_REAL * e1_force_RK(i,j,ispec,i_sls,3) + e1_force_RK(i,j,ispec,i_sls,4))
              endif

              e11_force_RK(i,j,ispec,i_sls,i_stage) = deltat * ((dux_dxl_n(i,j,ispec)-theta_n_u/TWO) * phinu2 - &
                                                                 e11(i,j,ispec,i_sls) * tauinvnu2)

              if( i_stage==1 .or. i_stage==2 .or. i_stage==3 ) then
                if( i_stage == 1) weight_rk = 0.5_CUSTOM_REAL
                if( i_stage == 2) weight_rk = 0.5_CUSTOM_REAL
                if( i_stage == 3) weight_rk = 1._CUSTOM_REAL

                if( i_stage==1) e11_initial_rk(i,j,ispec,i_sls) = e11(i,j,ispec,i_sls)
                e11(i,j,ispec,i_sls) = e11_initial_rk(i,j,ispec,i_sls) + weight_rk * e11_force_RK(i,j,ispec,i_sls,i_stage)
              else if( i_stage==4 ) then
                e11(i,j,ispec,i_sls) = e11_initial_rk(i,j,ispec,i_sls) + 1._CUSTOM_REAL / 6._CUSTOM_REAL * &
                                       (e11_force_RK(i,j,ispec,i_sls,1) + 2._CUSTOM_REAL * e11_force_RK(i,j,ispec,i_sls,2) + &
                                        2._CUSTOM_REAL * e11_force_RK(i,j,ispec,i_sls,3) + e11_force_RK(i,j,ispec,i_sls,4))
              endif

              e13_force_RK(i,j,ispec,i_sls,i_stage) = deltat * ((dux_dzl_n(i,j,ispec) + duz_dxl_n(i,j,ispec))*phinu2 - &
                                                                 e13(i,j,ispec,i_sls) * tauinvnu2)
              if( i_stage==1 .or. i_stage==2 .or. i_stage==3 ) then
                if( i_stage == 1) weight_rk = 0.5_CUSTOM_REAL
                if( i_stage == 2) weight_rk = 0.5_CUSTOM_REAL
                if( i_stage == 3) weight_rk = 1._CUSTOM_REAL

                if( i_stage==1) e13_initial_rk(i,j,ispec,i_sls) = e13(i,j,ispec,i_sls)
                e13(i,j,ispec,i_sls) = e13_initial_rk(i,j,ispec,i_sls) + weight_rk * e13_force_RK(i,j,ispec,i_sls,i_stage)
              else if( i_stage==4 ) then
                e13(i,j,ispec,i_sls) = e13_initial_rk(i,j,ispec,i_sls) + 1._CUSTOM_REAL / 6._CUSTOM_REAL * &
                                       (e13_force_RK(i,j,ispec,i_sls,1) + 2._CUSTOM_REAL * e13_force_RK(i,j,ispec,i_sls,2) + &
                                        2._CUSTOM_REAL * e13_force_RK(i,j,ispec,i_sls,3) + e13_force_RK(i,j,ispec,i_sls,4))
              endif
            endif
          enddo
        enddo
        enddo
      endif
    enddo
  endif
!!!! end of update memory variable in viscoelastic simulation

! this to avoid a warning at execution time about an undefined variable being used
! for the SH component in the case of a P-SV calculation, and vice versa
  sigma_xx = 0._CUSTOM_REAL; sigma_xy = 0._CUSTOM_REAL; sigma_xz = 0._CUSTOM_REAL
  sigma_zy = 0._CUSTOM_REAL; sigma_zz = 0._CUSTOM_REAL; sigma_zx = 0._CUSTOM_REAL

  if(  PML_BOUNDARY_CONDITIONS ) then
    accel_elastic_PML = 0._CUSTOM_REAL

    PML_dux_dxl = 0._CUSTOM_REAL; PML_dux_dzl = 0._CUSTOM_REAL
    PML_duz_dxl = 0._CUSTOM_REAL; PML_duz_dzl = 0._CUSTOM_REAL

    PML_dux_dxl_old = 0._CUSTOM_REAL; PML_dux_dzl_old = 0._CUSTOM_REAL
    PML_duz_dxl_old = 0._CUSTOM_REAL; PML_duz_dzl_old = 0._CUSTOM_REAL
  endif

  ifirstelem = 1
  ilastelem = nspec

  if( stage_time_scheme == 1 ) then
    time_n = (it-1) * deltat
    time_nsub1 = (it-2) * deltat
  else if( stage_time_scheme == 6 ) then
    time_n = (it-1) * deltat + c_LDDRK(i_stage) * deltat
  endif

  ! loop over spectral elements
  do ispec = ifirstelem,ilastelem
    tempx1(:,:) = 0._CUSTOM_REAL; tempy1(:,:) = 0._CUSTOM_REAL; tempz1(:,:) = 0._CUSTOM_REAL
    tempx2(:,:) = 0._CUSTOM_REAL; tempy2(:,:) = 0._CUSTOM_REAL; tempz2(:,:) = 0._CUSTOM_REAL
    tempx3(:,:) = 0._CUSTOM_REAL
    sigma_thetatheta(:,:) = 0._CUSTOM_REAL

    !--- elastic spectral element
    if( elastic(ispec) ) then
      ! get unrelaxed elastic parameters of current spectral element
      lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
      mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
      lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec))
      lambdalplusmul_unrelaxed_elastic = lambdal_unrelaxed_elastic + mul_unrelaxed_elastic

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX
          !--- if external medium, get elastic parameters of current grid point
          if( assign_external_model ) then
            cpl = vpext(i,j,ispec)
            csl = vsext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
            mul_unrelaxed_elastic = rhol*csl*csl
            lambdal_unrelaxed_elastic = rhol*cpl*cpl - TWO*mul_unrelaxed_elastic
            lambdaplus2mu_unrelaxed_elastic = lambdal_unrelaxed_elastic + TWO*mul_unrelaxed_elastic
            lambdalplusmul_unrelaxed_elastic = lambdal_unrelaxed_elastic + mul_unrelaxed_elastic
          endif

          ! derivative along x and along z
          dux_dxi = 0._CUSTOM_REAL;    duy_dxi = 0._CUSTOM_REAL;    duz_dxi = 0._CUSTOM_REAL
          dux_dgamma = 0._CUSTOM_REAL; duy_dgamma = 0._CUSTOM_REAL; duz_dgamma = 0._CUSTOM_REAL

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
            if( AXISYM ) then
              if (is_on_the_axis(ispec) ) then
                do k = 1,NGLJ
                  dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                  duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                  dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
                  duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i,k,ispec))*hprime_zz(j,k)
                enddo
              else
                do k = 1,NGLJ
                  dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
                  duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j,ispec))*hprime_xx(i,k)
                  dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
                  duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i,k,ispec))*hprime_zz(j,k)
                enddo
              endif
            else
              do k = 1,NGLLX
                dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
                duy_dxi = duy_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
                duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j,ispec))*hprime_xx(i,k)
                dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
                duy_dgamma = duy_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
                duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i,k,ispec))*hprime_zz(j,k)
              enddo
            endif


          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

          ! derivatives of displacement
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          duy_dxl = duy_dxi*xixl + duy_dgamma*gammaxl
          duy_dzl = duy_dxi*xizl + duy_dgamma*gammazl

          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

          if (AXISYM) then
            if(is_on_the_axis(ispec) .and. i == 1) then ! d_uz/dr=0 on the axis
              duz_dxl = 0.d0
            endif
          endif

          if( PML_BOUNDARY_CONDITIONS .and. is_PML(ispec) .and. nspec_PML > 0 ) then
            ispec_PML = spec_to_PML(ispec)
            CPML_region_local = region_CPML(ispec)
            kappa_x = K_x_store(i,j,ispec_PML)
            kappa_z = K_z_store(i,j,ispec_PML)
            d_x = d_x_store(i,j,ispec_PML)
            d_z = d_z_store(i,j,ispec_PML)
            alpha_x = alpha_x_store(i,j,ispec_PML)
            alpha_z = alpha_z_store(i,j,ispec_PML)
            beta_x = alpha_x + d_x / kappa_x
            beta_z = alpha_z + d_z / kappa_z

            PML_dux_dxl(i,j) = dux_dxl
            PML_dux_dzl(i,j) = dux_dzl
            PML_duz_dzl(i,j) = duz_dzl
            PML_duz_dxl(i,j) = duz_dxl

            ! derivative along x and along z
            dux_dxi_old = 0._CUSTOM_REAL;    duz_dxi_old = 0._CUSTOM_REAL
            dux_dgamma_old = 0._CUSTOM_REAL; duz_dgamma_old = 0._CUSTOM_REAL

            ! first double loop over GLL points to compute and store gradients
            ! we can merge the two loops because NGLLX == NGLLZ

            if( AXISYM ) then
              if (is_on_the_axis(ispec) ) then
                do k = 1,NGLJ
                  dux_dxi_old = dux_dxi_old + displ_elastic_old(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                  duz_dxi_old = duz_dxi_old + displ_elastic_old(3,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                  dux_dgamma_old = dux_dgamma_old + displ_elastic_old(1,ibool(i,k,ispec))*hprime_zz(j,k)
                  duz_dgamma_old = duz_dgamma_old + displ_elastic_old(3,ibool(i,k,ispec))*hprime_zz(j,k)
                enddo
              else
                do k = 1,NGLJ
                  dux_dxi_old = dux_dxi_old + displ_elastic_old(1,ibool(k,j,ispec))*hprime_xx(i,k)
                  duz_dxi_old = duz_dxi_old + displ_elastic_old(3,ibool(k,j,ispec))*hprime_xx(i,k)
                  dux_dgamma_old = dux_dgamma_old + displ_elastic_old(1,ibool(i,k,ispec))*hprime_zz(j,k)
                  duz_dgamma_old = duz_dgamma_old + displ_elastic_old(3,ibool(i,k,ispec))*hprime_zz(j,k)
                enddo
              endif
            else
              do k = 1,NGLLX
                dux_dxi_old = dux_dxi_old + displ_elastic_old(1,ibool(k,j,ispec))*hprime_xx(i,k)
                duz_dxi_old = duz_dxi_old + displ_elastic_old(3,ibool(k,j,ispec))*hprime_xx(i,k)
                dux_dgamma_old = dux_dgamma_old + displ_elastic_old(1,ibool(i,k,ispec))*hprime_zz(j,k)
                duz_dgamma_old = duz_dgamma_old + displ_elastic_old(3,ibool(i,k,ispec))*hprime_zz(j,k)
              enddo
            endif

            ! derivatives of displacement
            PML_dux_dxl_old(i,j) = dux_dxi_old*xixl + dux_dgamma_old*gammaxl !dux_dxl_old
            PML_dux_dzl_old(i,j) = dux_dxi_old*xizl + dux_dgamma_old*gammazl !dux_dzl_old
            PML_duz_dxl_old(i,j) = duz_dxi_old*xixl + duz_dgamma_old*gammaxl !duz_dxl_old
            PML_duz_dzl_old(i,j) = duz_dxi_old*xizl + duz_dgamma_old*gammazl !duz_dzl_old

            if (AXISYM .and. is_on_the_axis(ispec) .and. i == 1) then ! d_uz/dr=0 on the axis
              PML_duz_dxl_old(i,j) = 0.d0
            endif

!------------------------------------------------------------------------------
!---------------------------- LEFT & RIGHT ------------------------------------
!------------------------------------------------------------------------------
            call lik_parameter_computation(time_n,deltat,kappa_z,beta_z,alpha_z,kappa_x,beta_x,alpha_x,&
                                           CPML_region_local,31,A5,A6,A7,singularity_type_zx,bb_zx_1,bb_zx_2,&
                                           coef0_zx_1,coef1_zx_1,coef2_zx_1,coef0_zx_2,coef1_zx_2,coef2_zx_2)

            call lik_parameter_computation(time_n,deltat,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z,&
                                           CPML_region_local,13,A8,A9,A10,singularity_type_xz,bb_xz_1,bb_xz_2,&
                                           coef0_xz_1,coef1_xz_1,coef2_xz_1,coef0_xz_2,coef1_xz_2,coef2_xz_2)
            if( stage_time_scheme == 1 ) then
              if( ROTATE_PML_ACTIVATE ) then
                rmemory_dux_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_dux_dx(i,j,ispec_PML,1) + &
                                                  coef1_zx_1 * PML_dux_dxl(i,j) + coef2_zx_1 * PML_dux_dxl_old(i,j)
                rmemory_dux_dz(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_dux_dz(i,j,ispec_PML,1) + &
                                                  coef1_zx_1 * PML_dux_dzl(i,j) + coef2_zx_1 * PML_dux_dzl_old(i,j)
                rmemory_duz_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_duz_dx(i,j,ispec_PML,1) + &
                                                  coef1_zx_1 * PML_duz_dxl(i,j) + coef2_zx_1 * PML_duz_dxl_old(i,j)
                rmemory_duz_dz(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_duz_dz(i,j,ispec_PML,1) + &
                                                  coef1_zx_1 * PML_duz_dzl(i,j) + coef2_zx_1 * PML_duz_dzl_old(i,j)
                rmemory_dux_dx_prime(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_dux_dx_prime(i,j,ispec_PML,1) + &
                                                        coef1_xz_1 * PML_dux_dxl(i,j) + coef2_xz_1 * PML_dux_dxl_old(i,j)
                rmemory_dux_dz_prime(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_dux_dz_prime(i,j,ispec_PML,1) + &
                                                        coef1_xz_1 * PML_dux_dzl(i,j) + coef2_xz_1 * PML_dux_dzl_old(i,j)
                rmemory_duz_dx_prime(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_duz_dx_prime(i,j,ispec_PML,1) + &
                                                        coef1_xz_1 * PML_duz_dxl(i,j) + coef2_xz_1 * PML_duz_dxl_old(i,j)
                rmemory_duz_dz_prime(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_duz_dz_prime(i,j,ispec_PML,1) + &
                                                        coef1_xz_1 * PML_duz_dzl(i,j) + coef2_xz_1 * PML_duz_dzl_old(i,j)
                if( singularity_type_zx == 0 ) then
                  rmemory_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_dux_dx(i,j,ispec_PML,2) + &
                                                    coef1_zx_2 * PML_dux_dxl(i,j) + coef2_zx_2 * PML_dux_dxl_old(i,j)
                  rmemory_dux_dz(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_dux_dz(i,j,ispec_PML,2) + &
                                                    coef1_zx_2 * PML_dux_dzl(i,j) + coef2_zx_2 * PML_dux_dzl_old(i,j)
                  rmemory_duz_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_duz_dx(i,j,ispec_PML,2) + &
                                                    coef1_zx_2 * PML_duz_dxl(i,j) + coef2_zx_2 * PML_duz_dxl_old(i,j)
                  rmemory_duz_dz(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_duz_dz(i,j,ispec_PML,2) + &
                                                    coef1_zx_2 * PML_duz_dzl(i,j) + coef2_zx_2 * PML_duz_dzl_old(i,j)
                else
                  rmemory_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_dux_dx(i,j,ispec_PML,2) + &
                                                    coef1_zx_2 * time_n * PML_dux_dxl(i,j) + &
                                                    coef2_zx_2 * time_nsub1 * PML_dux_dxl_old(i,j)
                  rmemory_dux_dz(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_dux_dz(i,j,ispec_PML,2) + &
                                                    coef1_zx_2 * time_n * PML_dux_dzl(i,j) + &
                                                    coef2_zx_2 * time_nsub1 * PML_dux_dzl_old(i,j)
                  rmemory_duz_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_duz_dx(i,j,ispec_PML,2) + &
                                                    coef1_zx_2 * time_n * PML_duz_dxl(i,j) + &
                                                    coef2_zx_2 * time_nsub1 * PML_duz_dxl_old(i,j)
                  rmemory_duz_dz(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_duz_dz(i,j,ispec_PML,2) + &
                                                    coef1_zx_2 * time_n * PML_duz_dzl(i,j) + &
                                                    coef2_zx_2 * time_nsub1 * PML_duz_dzl_old(i,j)
                endif

                if( singularity_type_xz == 0 ) then
                  rmemory_dux_dx_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_dux_dx_prime(i,j,ispec_PML,2) + &
                                                          coef1_xz_2 * PML_dux_dxl(i,j) + coef2_xz_2 * PML_dux_dxl_old(i,j)
                  rmemory_dux_dz_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_dux_dz_prime(i,j,ispec_PML,2) + &
                                                          coef1_xz_2 * PML_dux_dzl(i,j) + coef2_xz_2 * PML_dux_dzl_old(i,j)
                  rmemory_duz_dx_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_duz_dx_prime(i,j,ispec_PML,2) + &
                                                          coef1_xz_2 * PML_duz_dxl(i,j) + coef2_xz_2 * PML_duz_dxl_old(i,j)
                  rmemory_duz_dz_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_duz_dz_prime(i,j,ispec_PML,2) + &
                                                          coef1_xz_2 * PML_duz_dzl(i,j) + coef2_xz_2 * PML_duz_dzl_old(i,j)
                else
                  rmemory_dux_dx_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_dux_dx_prime(i,j,ispec_PML,2) + &
                                                          coef1_xz_2 * time_n * PML_dux_dxl(i,j) + &
                                                          coef2_xz_2 * time_nsub1 * PML_dux_dxl_old(i,j)
                  rmemory_dux_dz_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_dux_dz_prime(i,j,ispec_PML,2) + &
                                                          coef1_xz_2 * time_n * PML_dux_dzl(i,j) + &
                                                          coef2_xz_2 * time_nsub1 * PML_dux_dzl_old(i,j)
                  rmemory_duz_dx_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_duz_dx_prime(i,j,ispec_PML,2) + &
                                                          coef1_xz_2 * time_n * PML_duz_dxl(i,j) + &
                                                          coef2_xz_2 * time_nsub1 * PML_duz_dxl_old(i,j)
                  rmemory_duz_dz_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_duz_dz_prime(i,j,ispec_PML,2) + &
                                                          coef1_xz_2 * time_n * PML_duz_dzl(i,j) + &
                                                          coef2_xz_2 * time_nsub1 * PML_duz_dzl_old(i,j)
                endif

              else
                rmemory_dux_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_dux_dx(i,j,ispec_PML,1) + &
                                                  coef1_zx_1 * PML_dux_dxl(i,j) + coef2_zx_1 * PML_dux_dxl_old(i,j)
                rmemory_duz_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_duz_dx(i,j,ispec_PML,1) + &
                                                  coef1_zx_1 * PML_duz_dxl(i,j) + coef2_zx_1 * PML_duz_dxl_old(i,j)
                if( singularity_type_zx == 0 ) then
                  rmemory_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_dux_dx(i,j,ispec_PML,2) + &
                                                    coef1_zx_2 * PML_dux_dxl(i,j) + coef2_zx_2 * PML_dux_dxl_old(i,j)
                  rmemory_duz_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_duz_dx(i,j,ispec_PML,2) + &
                                                    coef1_zx_2 * PML_duz_dxl(i,j) + coef2_zx_2 * PML_duz_dxl_old(i,j)
                else
                  rmemory_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_dux_dx(i,j,ispec_PML,2) + &
                                                    coef1_zx_2 * time_n * PML_dux_dxl(i,j) + &
                                                    coef2_zx_2 * time_nsub1 * PML_dux_dxl_old(i,j)
                  rmemory_duz_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_duz_dx(i,j,ispec_PML,2) + &
                                                    coef1_zx_2 * time_n * PML_duz_dxl(i,j) + &
                                                    coef2_zx_2 * time_nsub1 * PML_duz_dxl_old(i,j)
                endif

                rmemory_dux_dz(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_dux_dz(i,j,ispec_PML,1) + &
                                                  coef1_xz_1 * PML_dux_dzl(i,j) + coef2_xz_1 * PML_dux_dzl_old(i,j)
                rmemory_duz_dz(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_duz_dz(i,j,ispec_PML,1) + &
                                                  coef1_xz_1 * PML_duz_dzl(i,j) + coef2_xz_1 * PML_duz_dzl_old(i,j)
                if( singularity_type_xz == 0 ) then
                  rmemory_dux_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_dux_dz(i,j,ispec_PML,2) + &
                                                    coef1_xz_2 * PML_dux_dzl(i,j) + coef2_xz_2 * PML_dux_dzl_old(i,j)
                  rmemory_duz_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_duz_dz(i,j,ispec_PML,2) + &
                                                    coef1_xz_2 * PML_duz_dzl(i,j) + coef2_xz_2 * PML_duz_dzl_old(i,j)
                else
                  rmemory_dux_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_dux_dz(i,j,ispec_PML,2) + &
                                                    coef1_xz_2 * time_n * PML_dux_dzl(i,j) + &
                                                    coef2_xz_2 * time_nsub1 * PML_dux_dzl_old(i,j)
                  rmemory_duz_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_duz_dz(i,j,ispec_PML,2) + &
                                                    coef1_xz_2 * time_n * PML_duz_dzl(i,j) + &
                                                    coef2_xz_2 * time_nsub1 * PML_duz_dzl_old(i,j)
                endif

              endif
            endif

            if( stage_time_scheme == 6 ) then
              rmemory_dux_dx_LDDRK(i,j,ispec_PML,1) = alpha_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML,1) + &
                                                      deltat * (-bb_zx_1 * rmemory_dux_dx(i,j,ispec_PML,1) + PML_dux_dxl(i,j))
              rmemory_dux_dx(i,j,ispec_PML,1) = rmemory_dux_dx(i,j,ispec_PML,1) + &
                                                beta_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML,1)
              rmemory_duz_dx_LDDRK(i,j,ispec_PML,1) = alpha_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML,1) + &
                                                      deltat * (-bb_zx_1 * rmemory_duz_dx(i,j,ispec_PML,1) + PML_duz_dxl(i,j))
              rmemory_duz_dx(i,j,ispec_PML,1) = rmemory_duz_dx(i,j,ispec_PML,1) + &
                                                beta_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML,1)
              if( singularity_type_zx == 0 ) then
                rmemory_dux_dx_LDDRK(i,j,ispec_PML,2) = alpha_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML,2) + &
                                                        deltat * (-bb_zx_2 * rmemory_dux_dx(i,j,ispec_PML,2) + PML_dux_dxl(i,j))
                rmemory_dux_dx(i,j,ispec_PML,2) = rmemory_dux_dx(i,j,ispec_PML,2) + &
                                                  beta_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML,2)
                rmemory_duz_dx_LDDRK(i,j,ispec_PML,2) = alpha_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML,2) + &
                                                        deltat * (-bb_zx_2 * rmemory_duz_dx(i,j,ispec_PML,2) + PML_duz_dxl(i,j))
                rmemory_duz_dx(i,j,ispec_PML,2) = rmemory_duz_dx(i,j,ispec_PML,2) + &
                                                  beta_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML,2)
              else
                rmemory_dux_dx_LDDRK(i,j,ispec_PML,2) = alpha_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML,2) + &
                      deltat * (-bb_zx_2 * rmemory_dux_dx(i,j,ispec_PML,2) + PML_dux_dxl(i,j) * time_n)
                rmemory_dux_dx(i,j,ispec_PML,2) = rmemory_dux_dx(i,j,ispec_PML,2) + &
                                                  beta_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML,2)

                rmemory_duz_dx_LDDRK(i,j,ispec_PML,2) = alpha_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML,2) + &
                      deltat * (-bb_zx_2 * rmemory_duz_dx(i,j,ispec_PML,2) + PML_duz_dxl(i,j) * time_n)
                rmemory_duz_dx(i,j,ispec_PML,2) = rmemory_duz_dx(i,j,ispec_PML,2) + &
                                                  beta_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML,2)
              endif

              rmemory_dux_dz_LDDRK(i,j,ispec_PML,1) = alpha_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML,1) + &
                                                      deltat * (-bb_xz_1 * rmemory_dux_dz(i,j,ispec_PML,1) + PML_dux_dzl(i,j))
              rmemory_dux_dz(i,j,ispec_PML,1) = rmemory_dux_dz(i,j,ispec_PML,1) + &
                                                beta_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML,1)
              rmemory_duz_dz_LDDRK(i,j,ispec_PML,1) = alpha_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML,1) + &
                                                      deltat * (-bb_xz_1 * rmemory_duz_dz(i,j,ispec_PML,1) + PML_duz_dzl(i,j))
              rmemory_duz_dz(i,j,ispec_PML,1) = rmemory_duz_dz(i,j,ispec_PML,1) + &
                                                beta_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML,1)
              if( singularity_type_xz == 0 ) then
                rmemory_dux_dz_LDDRK(i,j,ispec_PML,2) = alpha_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML,2) + &
                                                        deltat * (-bb_xz_2 * rmemory_dux_dz(i,j,ispec_PML,2) + PML_dux_dzl(i,j))
                rmemory_dux_dz(i,j,ispec_PML,2) = rmemory_dux_dz(i,j,ispec_PML,2) + &
                                                  beta_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML,2)

                rmemory_duz_dz_LDDRK(i,j,ispec_PML,2) = alpha_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML,2) + &
                                                        deltat * (-bb_xz_2 * rmemory_duz_dz(i,j,ispec_PML,2) + PML_duz_dzl(i,j))
                rmemory_duz_dz(i,j,ispec_PML,2) = rmemory_duz_dz(i,j,ispec_PML,2) + &
                                                  beta_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML,2)
              else
                rmemory_dux_dz_LDDRK(i,j,ispec_PML,2) = alpha_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML,2) + &
                      deltat * (-bb_xz_2 * rmemory_dux_dz(i,j,ispec_PML,2) + PML_dux_dzl(i,j) * time_n)
                rmemory_dux_dz(i,j,ispec_PML,2) = rmemory_dux_dz(i,j,ispec_PML,2) + &
                                                  beta_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML,2)

                rmemory_duz_dz_LDDRK(i,j,ispec_PML,2) = alpha_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML,2) + &
                      deltat * (-bb_xz_2 * rmemory_duz_dz(i,j,ispec_PML,2) + PML_duz_dzl(i,j) * time_n)
                rmemory_duz_dz(i,j,ispec_PML,2) = rmemory_duz_dz(i,j,ispec_PML,2) + &
                                                  beta_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML,2)
              endif
            endif


            if( ROTATE_PML_ACTIVATE ) then
              dux_dxl = A5 * PML_dux_dxl(i,j)  + A6 * rmemory_dux_dx(i,j,ispec_PML,1) + A6 * rmemory_dux_dx(i,j,ispec_PML,2)
              dux_dzl = A5 * PML_dux_dzl(i,j)  + A6 * rmemory_dux_dz(i,j,ispec_PML,1) + A6 * rmemory_dux_dz(i,j,ispec_PML,2)
              duz_dxl = A5 * PML_duz_dxl(i,j)  + A6 * rmemory_duz_dx(i,j,ispec_PML,1) + A6 * rmemory_duz_dx(i,j,ispec_PML,2)
              duz_dzl = A5 * PML_duz_dzl(i,j)  + A6 * rmemory_duz_dz(i,j,ispec_PML,1) + A6 * rmemory_duz_dz(i,j,ispec_PML,2)

              dux_dxl_prime = A8 * PML_dux_dxl(i,j) + &
                              A9 * rmemory_dux_dx_prime(i,j,ispec_PML,1) + A10 * rmemory_dux_dx_prime(i,j,ispec_PML,2)
              dux_dzl_prime = A8 * PML_dux_dzl(i,j) + &
                              A9 * rmemory_dux_dz_prime(i,j,ispec_PML,1) + A10 * rmemory_dux_dz_prime(i,j,ispec_PML,2)
              duz_dxl_prime = A8 * PML_duz_dxl(i,j) + &
                              A9 * rmemory_duz_dx_prime(i,j,ispec_PML,1) + A10 * rmemory_duz_dx_prime(i,j,ispec_PML,2)
              duz_dzl_prime = A8 * PML_duz_dzl(i,j) + &
                              A9 * rmemory_duz_dz_prime(i,j,ispec_PML,1) + A10 * rmemory_duz_dz_prime(i,j,ispec_PML,2)
            else
              dux_dxl = A5 * PML_dux_dxl(i,j) + A6 * rmemory_dux_dx(i,j,ispec_PML,1) + A7 * rmemory_dux_dx(i,j,ispec_PML,2)
              duz_dxl = A5 * PML_duz_dxl(i,j) + A6 * rmemory_duz_dx(i,j,ispec_PML,1) + A7 * rmemory_duz_dx(i,j,ispec_PML,2)
              dux_dzl = A8 * PML_dux_dzl(i,j) + A9 * rmemory_dux_dz(i,j,ispec_PML,1) + A10 * rmemory_dux_dz(i,j,ispec_PML,2)
              duz_dzl = A8 * PML_duz_dzl(i,j) + A9 * rmemory_duz_dz(i,j,ispec_PML,1) + A10 * rmemory_duz_dz(i,j,ispec_PML,2)
            endif

            if (AXISYM .and. is_on_the_axis(ispec) .and. i == 1) then ! d_uz/dr=0 on the axis
              rmemory_duz_dx(i,j,ispec_PML,1) = 0.d0
              rmemory_duz_dx(i,j,ispec_PML,2) = 0.d0
              if (time_stepping_scheme /= 1) then
                rmemory_duz_dx_LDDRK(i,j,ispec_PML,1) = 0.d0
                rmemory_duz_dx_LDDRK(i,j,ispec_PML,2) = 0.d0
              endif
              if( ROTATE_PML_ACTIVATE ) then
                rmemory_duz_dx_prime(i,j,ispec_PML,1) = 0.d0
                rmemory_duz_dx_prime(i,j,ispec_PML,2) = 0.d0
              endif
            endif
          endif ! PML_BOUNDARY_CONDITIONS

          if (AXISYM .and. is_on_the_axis(ispec) .and. i == 1) then ! d_uz/dr=0 on the axis
            duz_dxl = 0.d0
            duz_dxl_prime = 0.d0
          endif

          ! compute stress tensor (include attenuation or anisotropy if needed)
          if( ATTENUATION_VISCOELASTIC_SOLID ) then
            ! attenuation is implemented following the memory variable formulation of
            ! J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
            ! vol. 58(1), p. 110-120 (1993). More details can be found in
            ! J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation in a linear
            ! viscoelastic medium, Geophysical Journal International, vol. 95, p. 597-611 (1988).

! When implementing viscoelasticity according to the Carcione 1993 paper, attenuation is
! non-causal rather than causal i.e. wave speed up instead of slowing down
! when attenuation is turned on. We fixed that issue (which is not incorrect but non traditional)
! by taking the unrelaxed state (infinite frequency) as a reference instead of the relaxed state (zero frequency)
! and also using equations in Carcione's 2007 book.
! See file doc/old_problem_attenuation_reference_Specfem2D_fixed_by_Xie_Zhinan.pdf
! and doc/how_we_modified_Carcione_1993_to_make_it_causal_and_include_the_missing_1_over_L_factor.pdf

! See also J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993) for two memory-variable mechanisms (page 112).

! and J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation
! in a linear viscoelastic medium, Geophysical Journal International,
! vol. 95, p. 597-611 (1988) for two memory-variable mechanisms (page 604).

            if( AXISYM ) then
              if (is_on_the_axis(ispec) ) then
                if (is_on_the_axis(ispec) .and. i == 1) then ! First GLJ point
                  sigma_xx = 0._CUSTOM_REAL
                  sigma_zz = 0._CUSTOM_REAL
                  sigma_thetatheta(i,j) = 0._CUSTOM_REAL
                  xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                  r_xiplus1(i,j) = xxi
                  do k = 1,NGLJ
                    sigma_xx = sigma_xx + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                    sigma_zz = sigma_zz + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                    sigma_thetatheta(i,j) = sigma_thetatheta(i,j) + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                  enddo
                  sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl &
                             + lambdal_unrelaxed_elastic*sigma_xx/xxi
                  sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                             + lambdal_unrelaxed_elastic*sigma_zz/xxi
                  sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
                  sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                     + lambdaplus2mu_unrelaxed_elastic*sigma_thetatheta(i,j)/xxi
                else ! Not first GLJ point
                  sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl &
                             + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                  sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                             + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                  sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
                  sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                          + lambdaplus2mu_unrelaxed_elastic &
                                          * displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                  r_xiplus1(i,j) = coord(1,ibool(i,j,ispec))/(xiglj(i)+ONE)
                endif
              else ! Not on the axis
                sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl &
                           + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                           + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
                sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                        + lambdaplus2mu_unrelaxed_elastic &
                                        * displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
              endif
            else ! Not axisym
              sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
              sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
              sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl
            endif

            ! add the memory variables using the relaxed parameters (Carcione 2007 page 125)
            ! beware: there is a bug in Carcione's equation (2c) of his 1993 paper for sigma_zz, we fixed it in the code below.

! When implementing viscoelasticity according to the Carcione 1993 paper, attenuation is
! non-causal rather than causal i.e. wave speed up instead of slowing down
! when attenuation is turned on. We fixed that issue (which is not incorrect but non traditional)
! by taking the unrelaxed state (infinite frequency) as a reference instead of the relaxed state (zero frequency)
! and also using equations in Carcione's 2007 book.
! See file doc/old_problem_attenuation_reference_Specfem2D_fixed_by_Xie_Zhinan.pdf
! and doc/how_we_modified_Carcione_1993_to_make_it_causal_and_include_the_missing_1_over_L_factor.pdf

! See also J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993) for two memory-variable mechanisms (page 112).

! and J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation
! in a linear viscoelastic medium, Geophysical Journal International,
! vol. 95, p. 597-611 (1988) for two memory-variable mechanisms (page 604).

            e1_sum = 0._CUSTOM_REAL; e11_sum = 0._CUSTOM_REAL;  e13_sum = 0._CUSTOM_REAL
            do i_sls = 1,N_SLS
              e1_sum = e1_sum + e1(i,j,ispec,i_sls)
              e11_sum = e11_sum + e11(i,j,ispec,i_sls)
              e13_sum = e13_sum + e13(i,j,ispec,i_sls)
            enddo

! use the right formula with 1/N included
! i.e. use the unrelaxed moduli here (see Carcione's book, third edition, equation (3.189))
            sigma_xx = sigma_xx + lambdalplusmul_unrelaxed_elastic * e1_sum + TWO * mul_unrelaxed_elastic * e11_sum
            sigma_xz = sigma_xz + mul_unrelaxed_elastic * e13_sum
            sigma_zz = sigma_zz + lambdalplusmul_unrelaxed_elastic * e1_sum - TWO * mul_unrelaxed_elastic * e11_sum
            sigma_zx = sigma_xz

            if( PML_BOUNDARY_CONDITIONS .and. is_PML(ispec) ) then
! PML currently has no support for viscoelasticity, use the elastic formula instead
              sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*PML_duz_dzl(i,j)
              sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*PML_dux_dxl(i,j)
              sigma_zx = mul_unrelaxed_elastic * (PML_duz_dxl(i,j) + dux_dzl)
              sigma_xz = mul_unrelaxed_elastic * (PML_dux_dzl(i,j) + duz_dxl)
            endif
          else
            ! no attenuation

            if( AXISYM ) then
              if (is_on_the_axis(ispec) ) then
                if (is_on_the_axis(ispec) .and. i == 1) then ! First GLJ point
                  sigma_xx = 0._CUSTOM_REAL
                  sigma_zz = 0._CUSTOM_REAL
                  sigma_thetatheta(i,j) = 0._CUSTOM_REAL
                  xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                  r_xiplus1(i,j) = xxi
                  do k = 1,NGLJ
                    sigma_xx = sigma_xx + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                    sigma_zz = sigma_zz + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                    sigma_thetatheta(i,j) = sigma_thetatheta(i,j) + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                  enddo
                  sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl &
                             + lambdal_unrelaxed_elastic*sigma_xx/xxi
                  sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                             + lambdal_unrelaxed_elastic*sigma_zz/xxi
                  sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
                  sigma_zx = sigma_xz
                  sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                     + lambdaplus2mu_unrelaxed_elastic*sigma_thetatheta(i,j)/xxi
                else ! Not first GLJ point
                  sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl &
                             + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                  sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                             + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                  sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
                  sigma_zx = sigma_xz
                  sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                          + lambdaplus2mu_unrelaxed_elastic &
                                          * displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                  r_xiplus1(i,j) = coord(1,ibool(i,j,ispec))/(xiglj(i)+ONE)
                endif
              else ! Not on the axis
                sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl &
                           + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                           + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
                sigma_zx = sigma_xz
                sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                        + lambdaplus2mu_unrelaxed_elastic &
                                        * displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
              endif
            else ! Not axisym
              sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
              sigma_xy = mul_unrelaxed_elastic*duy_dxl
              sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
              sigma_zy = mul_unrelaxed_elastic*duy_dzl
              sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl
              sigma_zx = sigma_xz
            endif

            if( PML_BOUNDARY_CONDITIONS .and. is_PML(ispec) .and. nspec_PML > 0 ) then
              if( ROTATE_PML_ACTIVATE ) then
                theta = -ROTATE_PML_ANGLE/180._CUSTOM_REAL*Pi
                if( it==1)write(*,*)theta,ROTATE_PML_ACTIVATE,cos(theta),sin(theta)
                ct=cos(theta)
                st=sin(theta)
                sigma_xx_prime = lambdaplus2mu_unrelaxed_elastic*(ct**2*dux_dxl+ct*st*duz_dxl+ct*st*dux_dzl+st**2*duz_dzl) + &
                                 lambdal_unrelaxed_elastic*(st**2*PML_dux_dxl(i,j) - ct*st*PML_duz_dxl(i,j) - &
                                                            ct*st*PML_dux_dzl(i,j) + ct**2*PML_duz_dzl(i,j))

                sigma_xz_prime = mul_unrelaxed_elastic * (-ct*st*dux_dxl+ct**2*duz_dxl-st**2*dux_dzl+ct*st*duz_dzl) + &
                                 mul_unrelaxed_elastic * (-ct*st*PML_dux_dxl(i,j) - st**2*PML_duz_dxl(i,j) + &
                                                          ct**2*PML_dux_dzl(i,j) + ct*st*PML_duz_dzl(i,j))

                sigma_zx_prime = mul_unrelaxed_elastic * (-ct*st*PML_dux_dxl(i,j) + ct**2*PML_duz_dxl(i,j) - &
                                                          st**2*PML_dux_dzl(i,j) + ct*st*PML_duz_dzl(i,j)) + &
                                 mul_unrelaxed_elastic * (-ct*st*dux_dxl_prime - st**2*duz_dxl_prime + &
                                                          ct**2*dux_dzl_prime + ct*st*duz_dzl_prime)

                sigma_zz_prime = lambdaplus2mu_unrelaxed_elastic*(st**2*dux_dxl_prime - ct*st*duz_dxl_prime - &
                                                                  ct*st*dux_dzl_prime + ct**2*duz_dzl_prime) + &
                                 lambdal_unrelaxed_elastic*(ct**2*PML_dux_dxl(i,j) + ct*st*PML_duz_dxl(i,j) + &
                                                            ct*st*PML_dux_dzl(i,j) + st**2*PML_duz_dzl(i,j))

                sigma_xx = ct**2*sigma_xx_prime-ct*st*sigma_xz_prime-ct*st*sigma_zx_prime+st**2*sigma_zz_prime
                sigma_xz = ct*st*sigma_xx_prime+ct**2*sigma_xz_prime-st**2*sigma_zx_prime-ct*st*sigma_zz_prime
                sigma_zx = ct*st*sigma_xx_prime-st**2*sigma_xz_prime+ct**2*sigma_zx_prime-ct*st*sigma_zz_prime
                sigma_zz = st**2*sigma_xx_prime+ct*st*sigma_xz_prime+ct*st*sigma_zx_prime+ct**2*sigma_zz_prime
              else
                sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*PML_duz_dzl(i,j)
                sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*PML_dux_dxl(i,j)
                sigma_zx = mul_unrelaxed_elastic * (PML_duz_dxl(i,j) + dux_dzl)
                sigma_xz = mul_unrelaxed_elastic * (PML_dux_dzl(i,j) + duz_dxl)
              endif
            endif
          endif

          ! full anisotropy
          if( anisotropic(ispec) ) then
            if( assign_external_model ) then
              c11 = c11ext(i,j,ispec)
              c13 = c13ext(i,j,ispec)
              c15 = c15ext(i,j,ispec)
              c33 = c33ext(i,j,ispec)
              c35 = c35ext(i,j,ispec)
              c55 = c55ext(i,j,ispec)
              c12 = c12ext(i,j,ispec)
              c23 = c23ext(i,j,ispec)
              c25 = c25ext(i,j,ispec)
              if ( AXISYM ) then
                c22 = c22ext(i,j,ispec) ! This variable is used for axisym simulations only
              endif
            else
              c11 = anisotropy(1,kmato(ispec))
              c13 = anisotropy(2,kmato(ispec))
              c15 = anisotropy(3,kmato(ispec))
              c33 = anisotropy(4,kmato(ispec))
              c35 = anisotropy(5,kmato(ispec))
              c55 = anisotropy(6,kmato(ispec))
              c12 = anisotropy(7,kmato(ispec))
              c23 = anisotropy(8,kmato(ispec))
              c25 = anisotropy(9,kmato(ispec))
              if ( AXISYM ) then
                c22 = anisotropy(10,kmato(ispec)) ! This variable is used for axisym simulations only
              endif
            endif

            ! implement anisotropy in 2D
            if ( AXISYM ) then
              if ( is_on_the_axis(ispec) ) then
                if ( i == 1 ) then ! first GLJ point, on the axis
                  sigma_xx = 0._CUSTOM_REAL
                  sigma_zz = 0._CUSTOM_REAL
                  sigma_thetatheta(i,j) = 0._CUSTOM_REAL
                  xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                  r_xiplus1(i,j) = xxi
                  do k = 1,NGLJ ! Compute the sum
                    sigma_xx = sigma_xx + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                    sigma_zz = sigma_zz + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                    sigma_thetatheta(i,j) = sigma_thetatheta(i,j) + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                  enddo
                  sigma_xx = c11*dux_dxl + c13*duz_dzl + c12*sigma_xx/xxi
                  sigma_zz = c13*dux_dxl + c33*duz_dzl + c23*sigma_zz/xxi
                  sigma_xz = c15*dux_dxl + c35*duz_dzl + c55*(duz_dxl + dux_dzl)
                  sigma_thetatheta(i,j) = c12*dux_dxl + c23*duz_dzl + c22*sigma_thetatheta(i,j)/xxi

                else ! first GLJ point but not on the axis
                  sigma_xx = c11*dux_dxl + c13*duz_dzl + c12*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                  sigma_zz = c13*dux_dxl + c33*duz_dzl + c23*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                  sigma_xz = c15*dux_dxl + c35*duz_dzl + c55*(duz_dxl + dux_dzl)
                  sigma_thetatheta(i,j) = c12*dux_dxl + c23*duz_dzl + &
                                          c22*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                endif
              else ! axisym but not on the axis
                sigma_xx = c11*dux_dxl + c13*duz_dzl + c12*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                sigma_zz = c13*dux_dxl + c33*duz_dzl + c23*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                sigma_xz = c15*dux_dxl + c35*duz_dzl + c55*(duz_dxl + dux_dzl)
                sigma_thetatheta(i,j) = c12*dux_dxl + c23*duz_dzl + c22*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
              endif
            else ! not AXISYM
              sigma_xx = c11*dux_dxl + c13*duz_dzl + c15*(duz_dxl + dux_dzl)
              sigma_zz = c13*dux_dxl + c33*duz_dzl + c35*(duz_dxl + dux_dzl)
              sigma_xz = c15*dux_dxl + c35*duz_dzl + c55*(duz_dxl + dux_dzl)
            endif
          endif

          ! weak formulation term based on stress tensor (non-symmetric form)
          ! also add GLL integration weights
          jacobianl = jacobian(i,j,ispec)

          !! AB AB with the notations of Komatitsch & Tromp 1999 (with 3 -> 2) :
          ! tempx1(i,j) = w.J.F_{11}^{ij}
          ! tempz1(i,j) = w.J.F_{21}^{ij}
          ! tempx2(i,j) = w.J.F_{12}^{ij}
          ! tempz2(i,j) = w.J.F_{22}^{ij}

          if (AXISYM) then
            if (is_on_the_axis(ispec) ) then
              ! This is normal, we always add a contribution depending on the value on the axis
              ! i.e. we purposely sum something at point (i,j) with something at point (1,j)
              tempx3(i,j) = wzgll(j)*jacobian(1,j,ispec)*sigma_thetatheta(1,j)*hprimeBarwglj_xx(1,i)

              ! not first GLJ point
              if ( abs(coord(1,ibool(i,j,ispec))) > TINYVAL ) then
                if ( i == 1 ) then
                  write(*,*) "Element number:",ispec
                  call exit_MPI("error: an axial element is rotated. The code should have been stopped before. Check that your &
                   &coordinates are >> TINYVAL. Maybe you should also have a look to &
                   &doc/problematic_case_that_we_exclude_for_axisymmetric.pdf")
                endif
                tempx3(i,j) = tempx3(i,j) + wzgll(j)*wxglj(i)*jacobian(i,j,ispec) &
                              * sigma_thetatheta(i,j)/(xiglj(i)+ONE) ! this goes to accel_x
              endif

              tempx2(i,j) = r_xiplus1(i,j)*wxglj(i)*jacobianl &
                            * (sigma_xx*gammaxl+sigma_zx*gammazl) ! this goes to accel_x
              tempz2(i,j) = r_xiplus1(i,j)*wxglj(i)*jacobianl &
                            * (sigma_xz*gammaxl+sigma_zz*gammazl) ! this goes to accel_z
              tempx1(i,j) = r_xiplus1(i,j)*wzgll(j)*jacobianl &
                            * (sigma_xx*xixl+sigma_zx*xizl) ! this goes to accel_x
              tempz1(i,j) = r_xiplus1(i,j)*wzgll(j)*jacobianl &
                            * (sigma_xz*xixl+sigma_zz*xizl) ! this goes to accel_z
            else ! axisym but not on the axis
              tempx2(i,j) = coord(1,ibool(i,j,ispec))*wxgll(i)*jacobianl &
                            *(sigma_xx*gammaxl+sigma_zx*gammazl) ! this goes to accel_x
              tempz2(i,j) = coord(1,ibool(i,j,ispec))*wxgll(i)*jacobianl &
                            *(sigma_xz*gammaxl+sigma_zz*gammazl) ! this goes to accel_z
              tempx1(i,j) = coord(1,ibool(i,j,ispec))*wzgll(j)*jacobianl &
                            *(sigma_xx*xixl+sigma_zx*xizl) ! this goes to accel_x
              tempz1(i,j) = coord(1,ibool(i,j,ispec))*wzgll(j)*jacobianl &
                            *(sigma_xz*xixl+sigma_zz*xizl) ! this goes to accel_z
              tempx3(i,j) = wxgll(i)*wzgll(j)*jacobianl*sigma_thetatheta(i,j) ! this goes to accel_x
            endif
          else

            tempx1(i,j) = wzgll(j)*jacobianl*(sigma_xx*xixl+sigma_zx*xizl) ! this goes to accel_x
            tempy1(i,j) = wzgll(j)*jacobianl*(sigma_xy*xixl+sigma_zy*xizl) ! this goes to accel_y
            tempz1(i,j) = wzgll(j)*jacobianl*(sigma_xz*xixl+sigma_zz*xizl) ! this goes to accel_z

            tempx2(i,j) = wxgll(i)*jacobianl*(sigma_xx*gammaxl+sigma_zx*gammazl) ! this goes to accel_x
            tempy2(i,j) = wxgll(i)*jacobianl*(sigma_xy*gammaxl+sigma_zy*gammazl) ! this goes to accel_y
            tempz2(i,j) = wxgll(i)*jacobianl*(sigma_xz*gammaxl+sigma_zz*gammazl) ! this goes to accel_z

          endif
        enddo
      enddo  ! end of the loops on the collocation points i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! update the displacement memory variable
      if( is_PML(ispec) .and. PML_BOUNDARY_CONDITIONS .and. nspec_PML > 0 ) then
        ispec_PML=spec_to_PML(ispec)
        CPML_region_local = region_CPML(ispec)
        do j = 1,NGLLZ
          do i = 1,NGLLX
            if(  assign_external_model ) then
              rhol = rhoext(i,j,ispec)
            else
              rhol = density(1,kmato(ispec))
            endif
            iglob=ibool(i,j,ispec)
            kappa_x = K_x_store(i,j,ispec_PML)
            kappa_z = K_z_store(i,j,ispec_PML)
            d_x = d_x_store(i,j,ispec_PML)
            d_z = d_z_store(i,j,ispec_PML)
            alpha_x = alpha_x_store(i,j,ispec_PML)
            alpha_z = alpha_z_store(i,j,ispec_PML)
            beta_x = alpha_x + d_x / kappa_x
            beta_z = alpha_z + d_z / kappa_z
            call l_parameter_computation(time_n,deltat,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                         CPML_region_local,A0,A1,A2,A3,A4,singularity_type,&
                                         bb_1,coef0_1,coef1_1,coef2_1,bb_2,coef0_2,coef1_2,coef2_2)

            if( stage_time_scheme == 1 ) then
              rmemory_displ_elastic(1,1,i,j,ispec_PML) = coef0_1 * rmemory_displ_elastic(1,1,i,j,ispec_PML) + &
                                                         coef1_1 * displ_elastic(1,iglob) + coef2_1 * displ_elastic_old(1,iglob)
              rmemory_displ_elastic(1,3,i,j,ispec_PML) = coef0_1 * rmemory_displ_elastic(1,3,i,j,ispec_PML) + &
                                                         coef1_1 * displ_elastic(3,iglob) + coef2_1 * displ_elastic_old(3,iglob)

              if( singularity_type == 0 ) then
                rmemory_displ_elastic(2,1,i,j,ispec_PML) = coef0_2 * rmemory_displ_elastic(2,1,i,j,ispec_PML) + &
                                                           coef1_2 * displ_elastic(1,iglob) + coef2_2 * displ_elastic_old(1,iglob)
                rmemory_displ_elastic(2,3,i,j,ispec_PML) = coef0_2 * rmemory_displ_elastic(2,3,i,j,ispec_PML) + &
                                                           coef1_2 * displ_elastic(3,iglob) + coef2_2 * displ_elastic_old(3,iglob)
              else
                rmemory_displ_elastic(2,1,i,j,ispec_PML) = coef0_2 * rmemory_displ_elastic(2,1,i,j,ispec_PML) + &
                                                           coef1_2 * time_n * displ_elastic(1,iglob) + &
                                                           coef2_2 * time_nsub1 * displ_elastic_old(1,iglob)
                rmemory_displ_elastic(2,3,i,j,ispec_PML) = coef0_2 * rmemory_displ_elastic(2,3,i,j,ispec_PML) + &
                                                           coef1_2 * time_n * displ_elastic(3,iglob) + &
                                                           coef2_2 * time_nsub1 * displ_elastic_old(3,iglob)
              endif
            endif

            if( stage_time_scheme == 6 ) then

              rmemory_displ_elastic_LDDRK(1,1,i,j,ispec_PML) = &
                    alpha_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(1,1,i,j,ispec_PML) + &
                    deltat * (-bb_1 * rmemory_displ_elastic(1,1,i,j,ispec_PML) + displ_elastic(1,iglob))
              rmemory_displ_elastic(1,1,i,j,ispec_PML) = rmemory_displ_elastic(1,1,i,j,ispec_PML) + &
                    beta_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(1,1,i,j,ispec_PML)

              rmemory_displ_elastic_LDDRK(1,3,i,j,ispec_PML) = &
                    alpha_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(1,3,i,j,ispec_PML) + &
                    deltat * (-bb_1 * rmemory_displ_elastic(1,3,i,j,ispec_PML) + displ_elastic(3,iglob))
              rmemory_displ_elastic(1,3,i,j,ispec_PML) = rmemory_displ_elastic(1,3,i,j,ispec_PML) + &
                    beta_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(1,3,i,j,ispec_PML)

              if( singularity_type == 0 ) then
                rmemory_displ_elastic_LDDRK(2,1,i,j,ispec_PML) = &
                      alpha_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,1,i,j,ispec_PML) + &
                      deltat * (-bb_2 * rmemory_displ_elastic(2,1,i,j,ispec_PML) + displ_elastic(1,iglob))
                rmemory_displ_elastic(2,1,i,j,ispec_PML) = rmemory_displ_elastic(2,1,i,j,ispec_PML) + &
                      beta_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,1,i,j,ispec_PML)

                rmemory_displ_elastic_LDDRK(2,3,i,j,ispec_PML) = &
                      alpha_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,3,i,j,ispec_PML) + &
                      deltat * (-bb_2 * rmemory_displ_elastic(2,3,i,j,ispec_PML) + displ_elastic(3,iglob))
                rmemory_displ_elastic(2,3,i,j,ispec_PML) = rmemory_displ_elastic(2,3,i,j,ispec_PML) + &
                      beta_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,3,i,j,ispec_PML)
              else
                rmemory_displ_elastic_LDDRK(2,1,i,j,ispec_PML) = &
                      alpha_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,1,i,j,ispec_PML) + &
                      deltat * (-bb_2 * rmemory_displ_elastic(2,1,i,j,ispec_PML) + displ_elastic(1,iglob) * time_n)
                rmemory_displ_elastic(2,1,i,j,ispec_PML) = rmemory_displ_elastic(2,1,i,j,ispec_PML) + &
                      beta_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,1,i,j,ispec_PML)

                rmemory_displ_elastic_LDDRK(2,3,i,j,ispec_PML) = &
                      alpha_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,3,i,j,ispec_PML) + &
                      deltat * (-bb_2 * rmemory_displ_elastic(2,3,i,j,ispec_PML) + displ_elastic(3,iglob) * time_n)
                rmemory_displ_elastic(2,3,i,j,ispec_PML) = rmemory_displ_elastic(2,3,i,j,ispec_PML) + &
                      beta_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,3,i,j,ispec_PML)
              endif

            endif

            if( AXISYM ) then
              if (is_on_the_axis(ispec) ) then
                accel_elastic_PML(1,i,j)= wxglj(i)*wzgll(j)*rhol*jacobian(i,j,ispec)*r_xiplus1(i,j) * &
                     ( A1 * veloc_elastic(1,iglob) + A2 * displ_elastic(1,iglob) + &
                       A3 * rmemory_displ_elastic(1,1,i,j,ispec_PML) + A4 * rmemory_displ_elastic(2,1,i,j,ispec_PML))
                accel_elastic_PML(3,i,j)= wxglj(i)*wzgll(j)*rhol*jacobian(i,j,ispec)*r_xiplus1(i,j) * &
                     ( A1 * veloc_elastic(3,iglob) + A2 * displ_elastic(3,iglob) + &
                       A3 * rmemory_displ_elastic(1,3,i,j,ispec_PML) + A4 * rmemory_displ_elastic(2,3,i,j,ispec_PML))
              else
                accel_elastic_PML(1,i,j)= wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec)*coord(1,ibool(i,j,ispec))* &
                     ( A1 * veloc_elastic(1,iglob) + A2 * displ_elastic(1,iglob) + &
                       A3 * rmemory_displ_elastic(1,1,i,j,ispec_PML) + A4 * rmemory_displ_elastic(2,1,i,j,ispec_PML))
                accel_elastic_PML(3,i,j)= wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec)*coord(1,ibool(i,j,ispec))* &
                     ( A1 * veloc_elastic(3,iglob) + A2 * displ_elastic(3,iglob) + &
                       A3 * rmemory_displ_elastic(1,3,i,j,ispec_PML) + A4 * rmemory_displ_elastic(2,3,i,j,ispec_PML))
              endif
            else

              accel_elastic_PML(1,i,j) = wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * &
                    ( A1 * veloc_elastic(1,iglob) + A2 * displ_elastic(1,iglob) + &
                      A3 * rmemory_displ_elastic(1,1,i,j,ispec_PML) + A4 * rmemory_displ_elastic(2,1,i,j,ispec_PML))
              accel_elastic_PML(3,i,j) = wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * &
                    ( A1 * veloc_elastic(3,iglob) + A2 * displ_elastic(3,iglob) + &
                      A3 * rmemory_displ_elastic(1,3,i,j,ispec_PML) + A4 * rmemory_displ_elastic(2,3,i,j,ispec_PML))
            endif
          enddo
        enddo
      endif ! update the displacement memory variable

      !
      ! second double-loop over GLL to compute all the terms
      !
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          ! along x direction and z direction
          ! and assemble the contributions
          ! we can merge the two loops because NGLLX == NGLLZ

          if (AXISYM ) then
            if (is_on_the_axis(ispec) ) then
              do k = 1,NGLJ
                accel_elastic(1,iglob) = accel_elastic(1,iglob) &
                                         - (tempx1(k,j)*hprimeBarwglj_xx(k,i) + tempx2(i,k)*hprimewgll_zz(k,j))
                accel_elastic(3,iglob) = accel_elastic(3,iglob) &
                                         - (tempz1(k,j)*hprimeBarwglj_xx(k,i) + tempz2(i,k)*hprimewgll_zz(k,j))
              enddo
              accel_elastic(1,iglob) = accel_elastic(1,iglob) - tempx3(i,j)
            else ! Axisym but not on the axis
              do k = 1,NGLLX
                accel_elastic(1,iglob) = accel_elastic(1,iglob) &
                                         - (tempx1(k,j)*hprimewgll_xx(k,i) + tempx2(i,k)*hprimewgll_zz(k,j))
                accel_elastic(3,iglob) = accel_elastic(3,iglob) &
                                         - (tempz1(k,j)*hprimewgll_xx(k,i) + tempz2(i,k)*hprimewgll_zz(k,j))
              enddo
              accel_elastic(1,iglob) = accel_elastic(1,iglob) - tempx3(i,j)
            endif
          else !if AXISYM == false
            do k = 1,NGLLX
              accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tempx1(k,j)*hprimewgll_xx(k,i) + tempx2(i,k)*hprimewgll_zz(k,j))
              accel_elastic(2,iglob) = accel_elastic(2,iglob) - (tempy1(k,j)*hprimewgll_xx(k,i) + tempy2(i,k)*hprimewgll_zz(k,j))
              accel_elastic(3,iglob) = accel_elastic(3,iglob) - (tempz1(k,j)*hprimewgll_xx(k,i) + tempz2(i,k)*hprimewgll_zz(k,j))
            enddo
          endif

          !!! PML_BOUNDARY_CONDITIONS
          if( is_PML(ispec) .and. PML_BOUNDARY_CONDITIONS ) then
            accel_elastic(1,iglob) = accel_elastic(1,iglob) - accel_elastic_PML(1,i,j)
            accel_elastic(3,iglob) = accel_elastic(3,iglob) - accel_elastic_PML(3,i,j)
          endif

        enddo
      enddo ! second loop over the GLL points
    endif ! end of test if elastic element
  enddo ! end of loop over all spectral elements

  !
  !--- Clayton-Engquist condition if elastic
  !
  if( STACEY_BOUNDARY_CONDITIONS ) then
    count_left=1
    count_right=1
    count_bottom=1

    do ispecabs = 1,nelemabs

      ispec = numabs(ispecabs)
      if( .not. elastic(ispec) ) cycle

      ! get elastic parameters of current spectral element
      lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
      mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
      lambdaplus2mu_unrelaxed_elastic = lambdal_unrelaxed_elastic + 2._CUSTOM_REAL * mul_unrelaxed_elastic
      rhol  = density(1,kmato(ispec))
      kappal  = lambdal_unrelaxed_elastic + TWO*mul_unrelaxed_elastic/3._CUSTOM_REAL
      cpl = sqrt((kappal + 4._CUSTOM_REAL*mul_unrelaxed_elastic/3._CUSTOM_REAL)/rhol)
      csl = sqrt(mul_unrelaxed_elastic/rhol)

      !--- left absorbing boundary
      if( codeabs(IEDGE4,ispecabs) ) then
        i = 1
        do j = 1,NGLLZ
          iglob = ibool(i,j,ispec)
          ! for analytical initial plane wave for Bielak's conditions
          ! left or right edge, horizontal normal vector
          if( add_Bielak_conditions .and. initialfield ) then
            if( .not.over_critical_angle ) then
              call compute_Bielak_conditions(coord,iglob,nglob,it,deltat,dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert, &
                          x0_source, z0_source, A_plane, B_plane, C_plane, anglesource(1), anglesource_refl, &
                          c_inc, c_refl, time_offset,f0)
              traction_x_t0 = lambdaplus2mu_unrelaxed_elastic * dxUx + lambdal_unrelaxed_elastic * dzUz
              traction_z_t0 = mul_unrelaxed_elastic * (dxUz + dzUx)
            else
              veloc_horiz=v0x_left(count_left)
              veloc_vert=v0z_left(count_left)
              traction_x_t0=t0x_left(count_left)
              traction_z_t0=t0z_left(count_left)
              count_left=count_left+1
            endif
          else
            veloc_horiz = 0._CUSTOM_REAL;   veloc_vert = 0._CUSTOM_REAL
            traction_x_t0 = 0._CUSTOM_REAL; traction_z_t0 = 0._CUSTOM_REAL
          endif

          ! external velocity model
          if( assign_external_model ) then
            cpl = vpext(i,j,ispec)
            csl = vsext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)
          nx = - zgamma / jacobian1D
          nz = + xgamma / jacobian1D

          weight = jacobian1D * wzgll(j)

          vx = veloc_elastic(1,iglob) - veloc_horiz
          vy = veloc_elastic(2,iglob)
          vz = veloc_elastic(3,iglob) - veloc_vert

          vn = nx*vx+nz*vz

          tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
          ty = rho_vs*vy
          tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

          displtx=0._CUSTOM_REAL; displtz=0._CUSTOM_REAL

          if( ADD_SPRING_TO_STACEY ) then
            displx = displ_elastic(1,iglob)
            disply = displ_elastic(2,iglob)
            displz = displ_elastic(3,iglob)

            spring_position=sqrt((coord(1,iglob)-x_center_spring)**2 + (coord(2,iglob)-z_center_spring)**2)

            displn = nx*displx+nz*displz

            displtx = lambdaplus2mu_unrelaxed_elastic / (2._CUSTOM_REAL*spring_position) * displn * nx + &
                      mul_unrelaxed_elastic / (2._CUSTOM_REAL*spring_position) * (displx-displn*nx)
            displty = mul_unrelaxed_elastic*disply
            displtz = lambdaplus2mu_unrelaxed_elastic / (2._CUSTOM_REAL*spring_position) * displn * nz + &
                      mul_unrelaxed_elastic / (2._CUSTOM_REAL*spring_position) * (displz-displn*nz)
          endif

          accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx + traction_x_t0+displtx)*weight
          accel_elastic(2,iglob) = accel_elastic(2,iglob) - ty*weight
          accel_elastic(3,iglob) = accel_elastic(3,iglob) - (tz + traction_z_t0+displtz)*weight

          if( SAVE_FORWARD .and. SIMULATION_TYPE ==1 ) then
            if( p_sv ) then !P-SV waves
              b_absorb_elastic_left(1,j,ib_left(ispecabs),it) = (tx + traction_x_t0)*weight
              b_absorb_elastic_left(3,j,ib_left(ispecabs),it) = (tz + traction_z_t0)*weight
            else !SH (membrane) waves
              b_absorb_elastic_left(2,j,ib_left(ispecabs),it) = ty*weight
            endif
          endif
        enddo
      endif  !  end of left absorbing boundary

      !--- right absorbing boundary
      if( codeabs(IEDGE2,ispecabs) ) then
        i = NGLLX
        do j = 1,NGLLZ
          ! Clayton-Engquist condition if elastic
          iglob = ibool(i,j,ispec)
          ! for analytical initial plane wave for Bielak's conditions
          ! left or right edge, horizontal normal vector
          if( add_Bielak_conditions .and. initialfield ) then
            if( .not.over_critical_angle ) then
              call compute_Bielak_conditions(coord,iglob,nglob,it,deltat,dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert, &
                          x0_source, z0_source, A_plane, B_plane, C_plane, anglesource(1), anglesource_refl, &
                          c_inc, c_refl, time_offset,f0)
              traction_x_t0 = lambdaplus2mu_unrelaxed_elastic * dxUx + lambdal_unrelaxed_elastic * dzUz
              traction_z_t0 = mul_unrelaxed_elastic*(dxUz + dzUx)
            else
              veloc_horiz=v0x_right(count_right)
              veloc_vert=v0z_right(count_right)
              traction_x_t0=t0x_right(count_right)
              traction_z_t0=t0z_right(count_right)
              count_right=count_right+1
            endif
          else
            veloc_horiz = 0._CUSTOM_REAL;   veloc_vert = 0._CUSTOM_REAL
            traction_x_t0 = 0._CUSTOM_REAL; traction_z_t0 = 0._CUSTOM_REAL
          endif

          ! external velocity model
          if( assign_external_model ) then
            cpl = vpext(i,j,ispec)
            csl = vsext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)
          nx = + zgamma / jacobian1D
          nz = - xgamma / jacobian1D

          weight = jacobian1D * wzgll(j)

          vx = veloc_elastic(1,iglob) - veloc_horiz
          vy = veloc_elastic(2,iglob)
          vz = veloc_elastic(3,iglob) - veloc_vert

          vn = nx*vx+nz*vz

          tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
          ty = rho_vs*vy
          tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

          displtx = 0._CUSTOM_REAL; displtz = 0._CUSTOM_REAL

          if( ADD_SPRING_TO_STACEY ) then
            displx = displ_elastic(1,iglob)
            disply = displ_elastic(2,iglob)
            displz = displ_elastic(3,iglob)

            spring_position=sqrt((coord(1,iglob)-x_center_spring)**2 + (coord(2,iglob)-z_center_spring)**2)

            displn = nx*displx+nz*displz

            displtx = lambdaplus2mu_unrelaxed_elastic / (2._CUSTOM_REAL*spring_position) * displn * nx + &
                       mul_unrelaxed_elastic / (2._CUSTOM_REAL*spring_position) * (displx-displn*nx)
            displty = mul_unrelaxed_elastic*disply
            displtz = lambdaplus2mu_unrelaxed_elastic / (2._CUSTOM_REAL*spring_position) * displn * nz + &
                      mul_unrelaxed_elastic / (2._CUSTOM_REAL*spring_position) * (displz-displn*nz)
          endif

          accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx - traction_x_t0+displtx)*weight
          accel_elastic(2,iglob) = accel_elastic(2,iglob) - ty*weight
          accel_elastic(3,iglob) = accel_elastic(3,iglob) - (tz - traction_z_t0+displtz)*weight

          if( SAVE_FORWARD .and. SIMULATION_TYPE ==1 ) then
            if( p_sv ) then !P-SV waves
              b_absorb_elastic_right(1,j,ib_right(ispecabs),it) = (tx - traction_x_t0)*weight
              b_absorb_elastic_right(3,j,ib_right(ispecabs),it) = (tz - traction_z_t0)*weight
            else! SH (membrane) waves
              b_absorb_elastic_right(2,j,ib_right(ispecabs),it) = ty*weight
            endif
          endif
        enddo
      endif  !  end of right absorbing boundary

      !--- bottom absorbing boundary
      if( codeabs(IEDGE1,ispecabs) ) then
        j = 1
!! DK DK not needed           ! exclude corners to make sure there is no contradiction on the normal
        ibegin = 1
        iend = NGLLX
!! DK DK not needed           if( codeabs(IEDGE4,ispecabs)) ibegin = 2
!! DK DK not needed           if( codeabs(IEDGE2,ispecabs)) iend = NGLLX-1
        do i = ibegin,iend
          ! Clayton-Engquist condition if elastic
          iglob = ibool(i,j,ispec)
          ! for analytical initial plane wave for Bielak's conditions
          ! top or bottom edge, vertical normal vector
          if( add_Bielak_conditions .and. initialfield ) then
            if( .not.over_critical_angle ) then
              call compute_Bielak_conditions(coord,iglob,nglob,it,deltat,dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert, &
                          x0_source, z0_source, A_plane, B_plane, C_plane, anglesource(1), anglesource_refl, &
                          c_inc, c_refl, time_offset,f0)
              traction_x_t0 = mul_unrelaxed_elastic * (dxUz + dzUx)
              traction_z_t0 = lambdal_unrelaxed_elastic * dxUx + lambdaplus2mu_unrelaxed_elastic * dzUz
            else
              veloc_horiz=v0x_bot(count_bottom)
              veloc_vert=v0z_bot(count_bottom)
              traction_x_t0=t0x_bot(count_bottom)
              traction_z_t0=t0z_bot(count_bottom)
              count_bottom=count_bottom+1
            endif
          else
            veloc_horiz = 0._CUSTOM_REAL;   veloc_vert = 0._CUSTOM_REAL
            traction_x_t0 = 0._CUSTOM_REAL; traction_z_t0 = 0._CUSTOM_REAL
          endif

          ! external velocity model
          if( assign_external_model ) then
            cpl = vpext(i,j,ispec)
            csl = vsext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)
          nx = + zxi / jacobian1D
          nz = - xxi / jacobian1D

          weight = jacobian1D * wxgll(i)

          vx = veloc_elastic(1,iglob) - veloc_horiz
          vy = veloc_elastic(2,iglob)
          vz = veloc_elastic(3,iglob) - veloc_vert

          vn = nx*vx+nz*vz

          tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
          ty = rho_vs*vy
          tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)
! exclude corners to make sure there is no contradiction on the normal
! for Stacey absorbing conditions but not for incident plane waves;
! thus subtract nothing i.e. zero in that case
          if( (codeabs_corner(1,ispecabs) .and. i == 1) .or. (codeabs_corner(2,ispecabs) .and. i == NGLLX) ) then
            tx = 0._CUSTOM_REAL; ty = 0._CUSTOM_REAL; tz = 0._CUSTOM_REAL
          endif

          displtx = 0._CUSTOM_REAL; displtz = 0._CUSTOM_REAL

          if( ADD_SPRING_TO_STACEY ) then
            displx = displ_elastic(1,iglob)
            disply = displ_elastic(2,iglob)
            displz = displ_elastic(3,iglob)

            spring_position=sqrt((coord(1,iglob)-x_center_spring)**2 + (coord(2,iglob)-z_center_spring)**2)

            displn = nx*displx+nz*displz

            displtx = lambdaplus2mu_unrelaxed_elastic / (2._CUSTOM_REAL*spring_position) * displn * nx + &
                      mul_unrelaxed_elastic / (2._CUSTOM_REAL*spring_position) * (displx-displn*nx)
            displty = mul_unrelaxed_elastic*disply
            displtz = lambdaplus2mu_unrelaxed_elastic / (2._CUSTOM_REAL*spring_position) * displn * nz + &
                      mul_unrelaxed_elastic / (2._CUSTOM_REAL*spring_position) * (displz-displn*nz)

            if( (codeabs(IEDGE4,ispecabs) .and. i == 1) .or. (codeabs(IEDGE2,ispecabs) .and. i == NGLLX) ) then
              displtx = 0._CUSTOM_REAL; displty = 0._CUSTOM_REAL; displtz = 0._CUSTOM_REAL
            endif
          endif

          accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx + traction_x_t0+displtx)*weight
          accel_elastic(2,iglob) = accel_elastic(2,iglob) - ty*weight
          accel_elastic(3,iglob) = accel_elastic(3,iglob) - (tz + traction_z_t0+displtz)*weight

          if( SAVE_FORWARD .and. SIMULATION_TYPE ==1 ) then
            if( p_sv ) then !P-SV waves
              b_absorb_elastic_bottom(1,i,ib_bottom(ispecabs),it) = (tx + traction_x_t0)*weight
              b_absorb_elastic_bottom(3,i,ib_bottom(ispecabs),it) = (tz + traction_z_t0)*weight
            else !SH (membrane) waves
              b_absorb_elastic_bottom(2,i,ib_bottom(ispecabs),it) = ty*weight
            endif
          endif
        enddo
      endif  !  end of bottom absorbing boundary

      !--- top absorbing boundary
      if( codeabs(IEDGE3,ispecabs) ) then
        j = NGLLZ
!! DK DK not needed           ! exclude corners to make sure there is no contradiction on the normal
        ibegin = 1
        iend = NGLLX
!! DK DK not needed           if( codeabs(IEDGE4,ispecabs)) ibegin = 2
!! DK DK not needed           if( codeabs(IEDGE2,ispecabs)) iend = NGLLX-1
        do i = ibegin,iend
          ! Clayton-Engquist condition if elastic
          iglob = ibool(i,j,ispec)
          ! for analytical initial plane wave for Bielak's conditions
          ! top or bottom edge, vertical normal vector
          if( add_Bielak_conditions .and. initialfield ) then
            call compute_Bielak_conditions(coord,iglob,nglob,it,deltat,dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert, &
                        x0_source, z0_source, A_plane, B_plane, C_plane, anglesource(1), anglesource_refl, &
                        c_inc, c_refl, time_offset,f0)
            traction_x_t0 = mul_unrelaxed_elastic * (dxUz + dzUx)
            traction_z_t0 = lambdal_unrelaxed_elastic * dxUx + lambdaplus2mu_unrelaxed_elastic * dzUz
          else
            veloc_horiz = 0._CUSTOM_REAL;   veloc_vert = 0._CUSTOM_REAL
            traction_x_t0 = 0._CUSTOM_REAL; traction_z_t0 = 0._CUSTOM_REAL
          endif

          ! external velocity model
          if( assign_external_model ) then
            cpl = vpext(i,j,ispec)
            csl = vsext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)
          nx = - zxi / jacobian1D
          nz = + xxi / jacobian1D

          weight = jacobian1D * wxgll(i)

          vx = veloc_elastic(1,iglob) - veloc_horiz
          vy = veloc_elastic(2,iglob)
          vz = veloc_elastic(3,iglob) - veloc_vert

          vn = nx*vx+nz*vz

          tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
          ty = rho_vs*vy
          tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)
! exclude corners to make sure there is no contradiction on the normal
! for Stacey absorbing conditions but not for incident plane waves;
! thus subtract nothing i.e. zero in that case
          if( (codeabs_corner(3,ispecabs) .and. i == 1) .or. (codeabs_corner(4,ispecabs) .and. i == NGLLX) ) then
            tx = 0._CUSTOM_REAL; ty = 0._CUSTOM_REAL; tz = 0._CUSTOM_REAL
          endif

          displtx = 0._CUSTOM_REAL; displtz = 0._CUSTOM_REAL

          if( ADD_SPRING_TO_STACEY ) then
            displx = displ_elastic(1,iglob)
            disply = displ_elastic(2,iglob)
            displz = displ_elastic(3,iglob)

            spring_position=sqrt((coord(1,iglob)-x_center_spring)**2 + (coord(2,iglob)-z_center_spring)**2)

            displn = nx*displx+nz*displz

            displtx = lambdaplus2mu_unrelaxed_elastic / (2._CUSTOM_REAL * spring_position) * displn * nx + &
                      mul_unrelaxed_elastic / (2._CUSTOM_REAL*spring_position) * (displx-displn*nx)
            displty = mul_unrelaxed_elastic * disply
            displtz = lambdaplus2mu_unrelaxed_elastic / (2._CUSTOM_REAL * spring_position) * displn * nz + &
                      mul_unrelaxed_elastic / (2._CUSTOM_REAL*spring_position) * (displz-displn*nz)

            if( (codeabs(IEDGE4,ispecabs) .and. i == 1) .or. (codeabs(IEDGE2,ispecabs) .and. i == NGLLX) ) then
              displtx = 0._CUSTOM_REAL; displty = 0._CUSTOM_REAL; displtz = 0._CUSTOM_REAL
            endif
          endif

          accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx + traction_x_t0+displtx)*weight
          accel_elastic(2,iglob) = accel_elastic(2,iglob) - ty*weight
          accel_elastic(3,iglob) = accel_elastic(3,iglob) - (tz + traction_z_t0+displtz)*weight

          if( SAVE_FORWARD .and. SIMULATION_TYPE ==1 ) then
            if( p_sv ) then !P-SV waves
              b_absorb_elastic_top(1,i,ib_top(ispecabs),it) = (tx- traction_x_t0)*weight
              b_absorb_elastic_top(3,i,ib_top(ispecabs),it) = (tz- traction_z_t0)*weight
            else !SH (membrane) waves
              b_absorb_elastic_top(2,i,ib_top(ispecabs),it) = ty*weight
            endif
          endif
        enddo
      endif  !  end of top absorbing boundary

    enddo
  endif  ! end of absorbing boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !set Dirichlet boundary condition on outer boundary of PML
  if(  PML_BOUNDARY_CONDITIONS .and. anyabs .and. nspec_PML > 0 ) then

! we have to put Dirichlet on the boundary of the PML
    do ispecabs = 1,nelemabs
      ispec = numabs(ispecabs)

      if (is_PML(ispec) ) then
        ispec_PML=spec_to_PML(ispec)

!--- left absorbing boundary
        if( codeabs(IEDGE4,ispecabs) ) then
          i = 1
          do j = 1,NGLLZ
            iglob = ibool(i,j,ispec)
            displ_elastic_old(:,iglob) = 0._CUSTOM_REAL
            displ_elastic(:,iglob) = 0._CUSTOM_REAL
            veloc_elastic(:,iglob) = 0._CUSTOM_REAL
            accel_elastic(:,iglob) = 0._CUSTOM_REAL
          enddo
        endif

!--- right absorbing boundary
        if( codeabs(IEDGE2,ispecabs) ) then
          i = NGLLX
          do j = 1,NGLLZ
            iglob = ibool(i,j,ispec)
            displ_elastic_old(:,iglob) = 0._CUSTOM_REAL
            displ_elastic(:,iglob) = 0._CUSTOM_REAL
            veloc_elastic(:,iglob) = 0._CUSTOM_REAL
            accel_elastic(:,iglob) = 0._CUSTOM_REAL
          enddo
        endif

!--- bottom absorbing boundary
        if( codeabs(IEDGE1,ispecabs) ) then
          j = 1
          ibegin = 1
          iend = NGLLX
          do i = ibegin,iend
            iglob = ibool(i,j,ispec)
            displ_elastic_old(:,iglob) = 0._CUSTOM_REAL
            displ_elastic(:,iglob) = 0._CUSTOM_REAL
            veloc_elastic(:,iglob) = 0._CUSTOM_REAL
            accel_elastic(:,iglob) = 0._CUSTOM_REAL
          enddo
        endif

!--- top absorbing boundary
        if( codeabs(IEDGE3,ispecabs) ) then
          j = NGLLZ
          ibegin = 1
          iend = NGLLX
          do i = ibegin,iend
            iglob = ibool(i,j,ispec)
            displ_elastic_old(:,iglob) = 0._CUSTOM_REAL
            displ_elastic(:,iglob) = 0._CUSTOM_REAL
            veloc_elastic(:,iglob) = 0._CUSTOM_REAL
            accel_elastic(:,iglob) = 0._CUSTOM_REAL
          enddo
        endif

      endif ! end of is_PML
    enddo ! end specabs loop
  endif  ! end of PML_BOUNDARY_CONDITIONS

end subroutine compute_forces_viscoelastic

