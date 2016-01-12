!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
! Copyright Universite de Pau et des Pays de l'Adour and CNRS, France.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
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

  subroutine compute_forces_poro_solid()

! compute forces for the solid poroelastic part

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ, &
    CPML_X_ONLY,CPML_Z_ONLY,IRIGHT,ILEFT,IBOTTOM,ITOP, &
    TWO,ONE,HALF,ZERO,FOUR_THIRDS, &
    IEDGE1,IEDGE2,IEDGE3,IEDGE4,ALPHA_LDDRK,BETA_LDDRK

  use specfem_par, only: nglob,nspec, &
                         ATTENUATION_VISCOELASTIC_SOLID,deltat, &
                         ibool,ispec_is_poroelastic, &
                         accels_poroelastic,displs_poroelastic,displw_poroelastic, &
                         displs_poroelastic_old, &
                         b_accels_poroelastic,b_displs_poroelastic,b_displw_poroelastic, &
                         xix,xiz,gammax,gammaz,jacobian, &
                         e11,e13,hprime_xx,hprimewgll_xx,hprime_zz,hprimewgll_zz,wxgll,wzgll, &
                         inv_tau_sigma_nu2,phi_nu2,Mu_nu2,N_SLS, &
                         mufr_k,B_k, &
                         SIMULATION_TYPE, &
                         e11_LDDRK,e13_LDDRK, &
                         e11_initial_rk,e13_initial_rk,e11_force_RK, e13_force_RK, &
                         stage_time_scheme,i_stage

  implicit none

  ! local variables
  integer :: ispec,i,j,k,iglob
  integer :: i_sls

  ! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: dwx_dxi,dwx_dgamma,dwz_dxi,dwz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  real(kind=CUSTOM_REAL) :: dwx_dxl,dwz_dxl,dwx_dzl,dwz_dzl
  real(kind=CUSTOM_REAL) :: b_dux_dxi,b_dux_dgamma,b_duz_dxi,b_duz_dgamma
  real(kind=CUSTOM_REAL) :: b_dux_dxl,b_duz_dxl,b_dux_dzl,b_duz_dzl
  real(kind=CUSTOM_REAL) :: dsxx,dsxz,dszz
  real(kind=CUSTOM_REAL) :: b_dsxx,b_dsxz,b_dszz
  real(kind=CUSTOM_REAL) :: b_dwx_dxi,b_dwx_dgamma,b_dwz_dxi,b_dwz_dgamma
  real(kind=CUSTOM_REAL) :: b_dwx_dxl,b_dwz_dxl,b_dwx_dzl,b_dwz_dzl
  real(kind=CUSTOM_REAL) :: dwxx,dwxz,dwzz
  real(kind=CUSTOM_REAL) :: b_dwxx,b_dwxz,b_dwzz
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_xz,sigma_zz
  real(kind=CUSTOM_REAL) :: sigmap
  real(kind=CUSTOM_REAL) :: b_sigma_xx,b_sigma_xz,b_sigma_zz
  real(kind=CUSTOM_REAL) :: b_sigmap

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2,tempz1,tempz2
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1p,tempx2p,tempz1p,tempz2p
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: b_tempx1,b_tempx2,b_tempz1,b_tempz2
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: b_tempx1p,b_tempx2p,b_tempz1p,b_tempz2p

! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

! material properties of the poroelastic medium
  real(kind=CUSTOM_REAL) :: mu_relaxed_viscoelastic,lambda_relaxed_viscoelastic,lambdalplus2mu_relaxed_viscoel
  real(kind=CUSTOM_REAL) :: mu_G,lambdal_G,lambdalplus2mul_G

  double precision :: phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar
  double precision :: D_biot,H_biot,C_biot,M_biot

! for attenuation
  real(kind=CUSTOM_REAL) :: phinu2,tauinvnu2,theta_n_u,theta_nsub1_u
  real(kind=CUSTOM_REAL) :: bb,coef0,coef1,coef2

  ! RK
  real(kind=CUSTOM_REAL) :: weight_rk

  real(kind=CUSTOM_REAL) :: e11_sum,e13_sum

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) ::dux_dxl_n,duz_dzl_n,duz_dxl_n,dux_dzl_n, &
                                                         !nsub1 denote discrete time step n-1
                                                         dux_dxl_nsub1,duz_dzl_nsub1,duz_dxl_nsub1,dux_dzl_nsub1

! implement attenuation
  if (ATTENUATION_VISCOELASTIC_SOLID) then

! compute Grad(displs_poroelastic) at time step n for attenuation
    call compute_gradient_attenuation(displs_poroelastic,dux_dxl_n,duz_dxl_n, &
           dux_dzl_n,duz_dzl_n,xix,xiz,gammax,gammaz,ibool,ispec_is_poroelastic,hprime_xx,hprime_zz,nspec,nglob)

! compute Grad(displs_poroelastic) at time step n-1 for attenuation
    call compute_gradient_attenuation(displs_poroelastic_old,dux_dxl_nsub1,duz_dxl_nsub1, &
           dux_dzl_nsub1,duz_dzl_nsub1,xix,xiz,gammax,gammaz,ibool,ispec_is_poroelastic,hprime_xx,hprime_zz,nspec,nglob)

! update memory variables with fourth-order Runge-Kutta time scheme for attenuation
! loop over spectral elements
  do ispec = 1,nspec

    if (ispec_is_poroelastic(ispec)) then
       do j = 1,NGLLZ; do i = 1,NGLLX
          theta_n_u = dux_dxl_n(i,j,ispec) + duz_dzl_n(i,j,ispec)
          theta_nsub1_u = dux_dxl_nsub1(i,j,ispec) + duz_dzl_nsub1(i,j,ispec)

          ! loop on all the standard linear solids
          do i_sls = 1,N_SLS
            phinu2 = phi_nu2(i,j,ispec,i_sls)
            tauinvnu2 = inv_tau_sigma_nu2(i,j,ispec,i_sls)

            ! update e1, e11, e13 in convolution formation with modified recursive convolution scheme on basis of
            ! second-order accurate convolution term calculation from equation (21) of
            ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
            ! Anisotropic-Medium PML for Vector FETD With Modified Basis Functions,
            ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)
            ! evolution e1 ! no need since we are just considering shear attenuation
            if (stage_time_scheme == 1) then
              bb = tauinvnu2; coef0 = exp(-bb * deltat)
              if (abs(bb) > 1e-5_CUSTOM_REAL) then
                 coef1 = (1._CUSTOM_REAL - exp(-bb * deltat / 2._CUSTOM_REAL)) / bb
                 coef2 = (1._CUSTOM_REAL - exp(-bb* deltat / 2._CUSTOM_REAL)) * exp(-bb * deltat / 2._CUSTOM_REAL)/ bb
              else
                 coef1 = deltat / 2._CUSTOM_REAL
                 coef2 = deltat / 2._CUSTOM_REAL
              endif

              e11(i,j,ispec,i_sls) = coef0 * e11(i,j,ispec,i_sls) + &
                                     phinu2 * (coef1 * (dux_dxl_n(i,j,ispec)-theta_n_u/TWO) + &
                                               coef2 * (dux_dxl_nsub1(i,j,ispec)-theta_nsub1_u/TWO))

              e13(i,j,ispec,i_sls) = coef0 * e13(i,j,ispec,i_sls) + &
                                     phinu2 * (coef1 * (dux_dzl_n(i,j,ispec) + duz_dxl_n(i,j,ispec)) + &
                                               coef2 * (dux_dzl_nsub1(i,j,ispec) + duz_dxl_nsub1(i,j,ispec)))
            endif

            ! update e1, e11, e13 in ADE formation with LDDRK scheme
            ! evolution e1 ! no need since we are just considering shear attenuation
            if (stage_time_scheme == 6) then
              e11_LDDRK(i,j,ispec,i_sls) = ALPHA_LDDRK(i_stage) * e11_LDDRK(i,j,ispec,i_sls) + &
                                           deltat * ((dux_dxl_n(i,j,ispec)-theta_n_u/TWO) * phinu2) - &
                                           deltat * (e11(i,j,ispec,i_sls) * tauinvnu2)
              e11(i,j,ispec,i_sls) = e11(i,j,ispec,i_sls)+BETA_LDDRK(i_stage)*e11_LDDRK(i,j,ispec,i_sls)

              e13_LDDRK(i,j,ispec,i_sls) = ALPHA_LDDRK(i_stage) * e13_LDDRK(i,j,ispec,i_sls) + &
                                           deltat * ((dux_dzl_n(i,j,ispec) + duz_dxl_n(i,j,ispec))*phinu2) - &
                                           deltat * (e13(i,j,ispec,i_sls) * tauinvnu2)
              e13(i,j,ispec,i_sls) = e13(i,j,ispec,i_sls)+BETA_LDDRK(i_stage) * e13_LDDRK(i,j,ispec,i_sls)
            endif

            ! update e1, e11, e13 in ADE formation with classical Runge-Kutta scheme
            ! evolution e1 ! no need since we are just considering shear attenuation
            if (stage_time_scheme == 4) then
              e11_force_RK(i,j,ispec,i_sls,i_stage) = deltat * ((dux_dxl_n(i,j,ispec)-theta_n_u/TWO) * phinu2 - &
                                                                 e11(i,j,ispec,i_sls) * tauinvnu2)
              if (i_stage==1 .or. i_stage==2 .or. i_stage==3) then
                if (i_stage == 1)weight_rk = 0.5_CUSTOM_REAL
                if (i_stage == 2)weight_rk = 0.5_CUSTOM_REAL
                if (i_stage == 3)weight_rk = 1._CUSTOM_REAL

                if (i_stage==1) e11_initial_rk(i,j,ispec,i_sls) = e11(i,j,ispec,i_sls)
                e11(i,j,ispec,i_sls) = e11_initial_rk(i,j,ispec,i_sls) + weight_rk * e11_force_RK(i,j,ispec,i_sls,i_stage)
              else if (i_stage==4) then
                e11(i,j,ispec,i_sls) = e11_initial_rk(i,j,ispec,i_sls) + 1._CUSTOM_REAL / 6._CUSTOM_REAL * &
                                       (e11_force_RK(i,j,ispec,i_sls,1) + 2._CUSTOM_REAL * e11_force_RK(i,j,ispec,i_sls,2) + &
                                        2._CUSTOM_REAL * e11_force_RK(i,j,ispec,i_sls,3) + e11_force_RK(i,j,ispec,i_sls,4))
              endif

              e13_force_RK(i,j,ispec,i_sls,i_stage) = deltat * ((dux_dzl_n(i,j,ispec) + duz_dxl_n(i,j,ispec))*phinu2 - &
                                                                 e13(i,j,ispec,i_sls) * tauinvnu2)
              if (i_stage==1 .or. i_stage==2 .or. i_stage==3) then
                if (i_stage == 1)weight_rk = 0.5_CUSTOM_REAL
                if (i_stage == 2)weight_rk = 0.5_CUSTOM_REAL
                if (i_stage == 3)weight_rk = 1._CUSTOM_REAL
                if (i_stage==1) e13_initial_rk(i,j,ispec,i_sls) = e13(i,j,ispec,i_sls)
                e13(i,j,ispec,i_sls) = e13_initial_rk(i,j,ispec,i_sls) + weight_rk * e13_force_RK(i,j,ispec,i_sls,i_stage)
              else if (i_stage==4) then
                e13(i,j,ispec,i_sls) = e13_initial_rk(i,j,ispec,i_sls) + 1._CUSTOM_REAL / 6._CUSTOM_REAL * &
                                       (e13_force_RK(i,j,ispec,i_sls,1) + 2._CUSTOM_REAL * e13_force_RK(i,j,ispec,i_sls,2) + &
                                        2._CUSTOM_REAL * e13_force_RK(i,j,ispec,i_sls,3) + e13_force_RK(i,j,ispec,i_sls,4))
              endif
            endif
          enddo
        enddo; enddo
     endif
   enddo

  endif ! end of test on attenuation

! loop over spectral elements
  do ispec = 1,nspec

!---
!--- poroelastic spectral element
!---

    if (ispec_is_poroelastic(ispec)) then

      ! get poroelastic parameters of current spectral element
      call get_poroelastic_material(ispec,phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar)

      ! Biot coefficients for the input phi
      call get_poroelastic_Biot_coeff(phi,kappa_s,kappa_f,kappa_fr,mu_fr,D_biot,H_biot,C_biot,M_biot)

      !The RHS has the form : div T -phi/c div T_f + phi/ceta_fk^-1.partial t w
      !where T = G:grad u_s + C_biot div w I
      !and T_f = C_biot div u_s I + M_biot div w I
      mu_G = mu_fr
      lambdal_G = H_biot - 2._CUSTOM_REAL*mu_fr
      lambdalplus2mul_G = lambdal_G + TWO*mu_G

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX
          ! derivative along x and along z for u_s and w
          dux_dxi = ZERO
          duz_dxi = ZERO

          dux_dgamma = ZERO
          duz_dgamma = ZERO

          dwx_dxi = ZERO
          dwz_dxi = ZERO

          dwx_dgamma = ZERO
          dwz_dgamma = ZERO

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displs_poroelastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displs_poroelastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displs_poroelastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displs_poroelastic(2,ibool(i,k,ispec))*hprime_zz(j,k)

            dwx_dxi = dwx_dxi + displw_poroelastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            dwz_dxi = dwz_dxi + displw_poroelastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            dwx_dgamma = dwx_dgamma + displw_poroelastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            dwz_dgamma = dwz_dgamma + displw_poroelastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
          enddo

          if (SIMULATION_TYPE == 3) then
            ! kernels calculation
            b_dux_dxi = ZERO
            b_duz_dxi = ZERO

            b_dux_dgamma = ZERO
            b_duz_dgamma = ZERO

            b_dwx_dxi = ZERO
            b_dwz_dxi = ZERO

            b_dwx_dgamma = ZERO
            b_dwz_dgamma = ZERO

            do k = 1,NGLLX
              b_dux_dxi = b_dux_dxi + b_displs_poroelastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
              b_duz_dxi = b_duz_dxi + b_displs_poroelastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
              b_dux_dgamma = b_dux_dgamma + b_displs_poroelastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
              b_duz_dgamma = b_duz_dgamma + b_displs_poroelastic(2,ibool(i,k,ispec))*hprime_zz(j,k)

              b_dwx_dxi = b_dwx_dxi + b_displw_poroelastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
              b_dwz_dxi = b_dwz_dxi + b_displw_poroelastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
              b_dwx_dgamma = b_dwx_dgamma + b_displw_poroelastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
              b_dwz_dgamma = b_dwz_dgamma + b_displw_poroelastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
            enddo
          endif

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

          ! derivatives of displacement
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

          dwx_dxl = dwx_dxi*xixl + dwx_dgamma*gammaxl
          dwx_dzl = dwx_dxi*xizl + dwx_dgamma*gammazl

          dwz_dxl = dwz_dxi*xixl + dwz_dgamma*gammaxl
          dwz_dzl = dwz_dxi*xizl + dwz_dgamma*gammazl

          if (SIMULATION_TYPE == 3) then ! kernels calculation
            b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
            b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl

            b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
            b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl

            b_dwx_dxl = b_dwx_dxi*xixl + b_dwx_dgamma*gammaxl
            b_dwx_dzl = b_dwx_dxi*xizl + b_dwx_dgamma*gammazl

            b_dwz_dxl = b_dwz_dxi*xixl + b_dwz_dgamma*gammaxl
            b_dwz_dzl = b_dwz_dxi*xizl + b_dwz_dgamma*gammazl
          endif

          ! compute stress tensor (include attenuation or anisotropy if needed)
          if (ATTENUATION_VISCOELASTIC_SOLID) then
! Dissipation only controlled by frame share attenuation in poroelastic (see Morency & Tromp, GJI 2008).
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

            ! compute unrelaxed elastic coefficients from formulas in Carcione 2007 page 125
            lambda_relaxed_viscoelastic = (lambdal_G + mu_G) - mu_G / Mu_nu2(i,j,ispec)
            mu_relaxed_viscoelastic = mu_G / Mu_nu2(i,j,ispec)
            lambdalplus2mu_relaxed_viscoel = lambda_relaxed_viscoelastic + TWO*mu_relaxed_viscoelastic

            ! compute the stress using the unrelaxed Lame parameters (Carcione 2007 page 125)
            sigma_xx = (lambdal_G + 2.0*mu_G)*dux_dxl + lambdal_G*duz_dzl + C_biot*(dwx_dxl + dwz_dzl)
            sigma_xz = mu_G*(duz_dxl + dux_dzl)
            sigma_zz = (lambdal_G + 2.0*mu_G)*duz_dzl + lambdal_G*dux_dxl + C_biot*(dwx_dxl + dwz_dzl)

            sigmap = C_biot*(dux_dxl + duz_dzl) + M_biot*(dwx_dxl + dwz_dzl)

            ! add the memory variables using the relaxed parameters (Carcione 2007 page 125)
            ! beware: there is a bug in Carcione's equation (2c) of his 1993 paper for sigma_zz, we fixed it in the code below
            e11_sum = 0._CUSTOM_REAL
            e13_sum = 0._CUSTOM_REAL

            do i_sls = 1,N_SLS
              e11_sum = e11_sum + e11(i,j,ispec,i_sls)
              e13_sum = e13_sum + e13(i,j,ispec,i_sls)
            enddo

            ! mu_G is the relaxed modulus. Note that it is defined as the
            ! frame modulus (in compute_forces_poro_solid.f90), which Christina Morency noted
            ! mu_fr, which is in her case equivalent to the solid phase shear
            ! modulus, and whose value is entered in Par_file for example
            sigma_xx = sigma_xx + TWO * mu_relaxed_viscoelastic * e11_sum
            sigma_xz = sigma_xz + mu_relaxed_viscoelastic * e13_sum
            sigma_zz = sigma_zz - TWO * mu_relaxed_viscoelastic * e11_sum

          else

            ! no attenuation
            sigma_xx = lambdalplus2mul_G*dux_dxl + lambdal_G*duz_dzl + C_biot*(dwx_dxl + dwz_dzl)
            sigma_xz = mu_G*(duz_dxl + dux_dzl)
            sigma_zz = lambdalplus2mul_G*duz_dzl + lambdal_G*dux_dxl + C_biot*(dwx_dxl + dwz_dzl)

            sigmap = C_biot*(dux_dxl + duz_dzl) + M_biot*(dwx_dxl + dwz_dzl)

            if (SIMULATION_TYPE == 3) then
              ! kernels calculation
              b_sigma_xx = lambdalplus2mul_G*b_dux_dxl + lambdal_G*b_duz_dzl + C_biot*(b_dwx_dxl + b_dwz_dzl)
              b_sigma_xz = mu_G*(b_duz_dxl + b_dux_dzl)
              b_sigma_zz = lambdalplus2mul_G*b_duz_dzl + lambdal_G*b_dux_dxl + C_biot*(b_dwx_dxl + b_dwz_dzl)

              b_sigmap = C_biot*(b_dux_dxl + b_duz_dzl) + M_biot*(b_dwx_dxl + b_dwz_dzl)
            endif
          endif

          ! kernels calculation
          if (SIMULATION_TYPE == 3) then
            iglob = ibool(i,j,ispec)
            dsxx =  dux_dxl
            dsxz = HALF * (duz_dxl + dux_dzl)
            dszz =  duz_dzl

            dwxx =  dwx_dxl
            dwxz = HALF * (dwz_dxl + dwx_dzl)
            dwzz =  dwz_dzl

            b_dsxx =  b_dux_dxl
            b_dsxz = HALF * (b_duz_dxl + b_dux_dzl)
            b_dszz =  b_duz_dzl

            b_dwxx =  b_dwx_dxl
            b_dwxz = HALF * (b_dwz_dxl + b_dwx_dzl)
            b_dwzz =  b_dwz_dzl

            B_k(iglob) = (dux_dxl + duz_dzl) *  (b_dux_dxl + b_duz_dzl) * (H_biot - FOUR_THIRDS * mu_fr)
            mufr_k(iglob) = (dsxx * b_dsxx + dszz * b_dszz + &
                            2._CUSTOM_REAL * dsxz * b_dsxz - &
                            1._CUSTOM_REAL/3._CUSTOM_REAL * (dux_dxl + duz_dzl) * (b_dux_dxl + b_duz_dzl) ) * mu_fr
          endif

          jacobianl = jacobian(i,j,ispec)

          ! weak formulation term based on stress tensor (symmetric form)
          ! also add GLL integration weights
          tempx1(i,j) = wzgll(j)*jacobianl*(sigma_xx*xixl+sigma_xz*xizl)
          tempz1(i,j) = wzgll(j)*jacobianl*(sigma_xz*xixl+sigma_zz*xizl)

          tempx2(i,j) = wxgll(i)*jacobianl*(sigma_xx*gammaxl+sigma_xz*gammazl)
          tempz2(i,j) = wxgll(i)*jacobianl*(sigma_xz*gammaxl+sigma_zz*gammazl)

          tempx1p(i,j) = wzgll(j)*jacobianl*sigmap*xixl
          tempz1p(i,j) = wzgll(j)*jacobianl*sigmap*xizl

          tempx2p(i,j) = wxgll(i)*jacobianl*sigmap*gammaxl
          tempz2p(i,j) = wxgll(i)*jacobianl*sigmap*gammazl

          if (SIMULATION_TYPE == 3) then ! kernels calculation
            b_tempx1(i,j) = wzgll(j)*jacobianl*(b_sigma_xx*xixl+b_sigma_xz*xizl)
            b_tempz1(i,j) = wzgll(j)*jacobianl*(b_sigma_xz*xixl+b_sigma_zz*xizl)

            b_tempx2(i,j) = wxgll(i)*jacobianl*(b_sigma_xx*gammaxl+b_sigma_xz*gammazl)
            b_tempz2(i,j) = wxgll(i)*jacobianl*(b_sigma_xz*gammaxl+b_sigma_zz*gammazl)

            b_tempx1p(i,j) = wzgll(j)*jacobianl*b_sigmap*xixl
            b_tempz1p(i,j) = wzgll(j)*jacobianl*b_sigmap*xizl

            b_tempx2p(i,j) = wxgll(i)*jacobianl*b_sigmap*gammaxl
            b_tempz2p(i,j) = wxgll(i)*jacobianl*b_sigmap*gammazl
          endif

        enddo
      enddo

!
! second double-loop over GLL to compute all the terms
!
      do j = 1,NGLLZ
        do i = 1,NGLLX

          iglob = ibool(i,j,ispec)

          ! along x direction and z direction
          ! and assemble the contributions
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) - ( (tempx1(k,j) - phi/tort*tempx1p(k,j)) &
                                          *hprimewgll_xx(k,i) + (tempx2(i,k) - phi/tort*tempx2p(i,k))*hprimewgll_zz(k,j))

            accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) - ( (tempz1(k,j) - phi/tort*tempz1p(k,j)) &
                                          *hprimewgll_xx(k,i) + (tempz2(i,k) - phi/tort*tempz2p(i,k))*hprimewgll_zz(k,j))
          enddo

          if (SIMULATION_TYPE == 3) then
            ! kernels calculation
            do k = 1,NGLLX
              b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) - ( (b_tempx1(k,j) - phi/tort*b_tempx1p(k,j)) &
                                              *hprimewgll_xx(k,i) + (b_tempx2(i,k) - phi/tort*b_tempx2p(i,k))*hprimewgll_zz(k,j))

              b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) - ( (b_tempz1(k,j) - phi/tort*b_tempz1p(k,j)) &
                                              *hprimewgll_xx(k,i) + (b_tempz2(i,k) - phi/tort*b_tempz2p(i,k))*hprimewgll_zz(k,j))
            enddo
          endif

        enddo ! second loop over the GLL points
      enddo

    endif ! end of test if poroelastic element

  enddo ! end of loop over all spectral elements

  end subroutine compute_forces_poro_solid

