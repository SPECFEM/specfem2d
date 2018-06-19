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

  subroutine compute_forces_viscoelastic_backward(b_accel_elastic,b_displ_elastic,b_displ_elastic_old, &
                                                  e1,e11,e13,iphase)

  ! compute forces for the elastic elements
  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,NDIM, &
    ONE,TWO,PI,TINYVAL,ALPHA_LDDRK,BETA_LDDRK

  use specfem_par, only: nglob,nspec,assign_external_model,P_SV, &
                         ATTENUATION_VISCOELASTIC,nspec_ATT,N_SLS, &
                         ibool,kmato,ispec_is_elastic, &
                         poroelastcoef,xix,xiz,gammax,gammaz, &
                         jacobian,vpext,vsext,rhoext,c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext, &
                         ispec_is_anisotropic,anisotropy, &
                         e1_LDDRK,e11_LDDRK,e13_LDDRK, &
                         e1_initial_rk,e11_initial_rk,e13_initial_rk,e1_force_RK, e11_force_RK, e13_force_RK, &
                         hprime_xx,hprimewgll_xx,hprime_zz,hprimewgll_zz,wxgll,wzgll, &
                         AXISYM,is_on_the_axis,hprimeBar_xx,hprimeBarwglj_xx,xiglj,wxglj, &
                         inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2,N_SLS, &
                         deltat,coord, &
                         time_stepping_scheme,i_stage,ispec_is_acoustic

  ! overlapping communication
  use specfem_par, only: nspec_inner_elastic,nspec_outer_elastic,phase_ispec_inner_elastic

  ! PML arrays
  use specfem_par, only: PML_BOUNDARY_CONDITIONS,ispec_is_PML
  ! CPML coefficients and memory variables
  ! for further optimization, we can exclude computation in PML region in the
  ! case of backward simulation

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob) :: b_accel_elastic,b_displ_elastic,b_displ_elastic_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_ATT,N_SLS) :: e1,e11,e13

  integer,intent(in) :: iphase

  !---
  !--- local variables
  !---

  integer :: ispec,i,j,k,iglob

  ! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_xy,sigma_xz,sigma_zy,sigma_zz,sigma_zx
  real(kind=CUSTOM_REAL) :: xxi

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2,tempz1,tempz2
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
    lambdalplusmul_unrelaxed_elastic,lambdaplus2mu_unrelaxed_elastic,cpl,csl,rhol

  ! for attenuation
  real(kind=CUSTOM_REAL) :: phinu1,phinu2,theta_n_u,theta_nsub1_u
  double precision :: tauinvnu1,tauinvnu2

  ! for anisotropy
  double precision ::  c11,c15,c13,c33,c35,c55,c12,c23,c25

  integer :: num_elements,ispec_p

  ! temp variable RK
  real(kind=CUSTOM_REAL) :: weight_rk

  !!!update memory variable in viscoelastic simulation
  if (iphase == 1) then

  if (ATTENUATION_VISCOELASTIC) then

    ! compute Grad(b_displ_elastic) at time step n for attenuation
    call compute_gradient_attenuation(b_displ_elastic,dux_dxl_n,duz_dxl_n, &
          dux_dzl_n,duz_dzl_n,xix,xiz,gammax,gammaz,ibool,ispec_is_elastic,hprime_xx,hprime_zz,nspec,nglob)

    ! compute Grad(disp_elastic_old) at time step n-1 for attenuation
    call compute_gradient_attenuation(b_displ_elastic_old,dux_dxl_nsub1,duz_dxl_nsub1, &
          dux_dzl_nsub1,duz_dzl_nsub1,xix,xiz,gammax,gammaz,ibool,ispec_is_elastic,hprime_xx,hprime_zz,nspec,nglob)

    ! loop over spectral elements
    do ispec = 1,nspec

      ! attenuation is not implemented in acoustic (i.e. fluid) media for now, only in viscoelastic (i.e. solid) media
      if (ispec_is_acoustic(ispec)) cycle

      if ((.not. PML_BOUNDARY_CONDITIONS) .or. (PML_BOUNDARY_CONDITIONS .and. (.not. ispec_is_PML(ispec)))) then
        do j = 1,NGLLZ
        do i = 1,NGLLX

          ! convention to indicate that Q = 9999 in that element i.e. that there is no viscoelasticity in that element
          if (inv_tau_sigma_nu1(i,j,ispec,1) < 0.) cycle

          theta_n_u = dux_dxl_n(i,j,ispec) + duz_dzl_n(i,j,ispec)
          theta_nsub1_u = dux_dxl_nsub1(i,j,ispec) + duz_dzl_nsub1(i,j,ispec)

          ! loop on all the standard linear solids
          do i_sls = 1,N_SLS
            phinu1 = phi_nu1(i,j,ispec,i_sls)
            tauinvnu1 = inv_tau_sigma_nu1(i,j,ispec,i_sls)
            phinu2 = phi_nu2(i,j,ispec,i_sls)
            tauinvnu2 = inv_tau_sigma_nu2(i,j,ispec,i_sls)

            ! update e1, e11, e13 in convolution formulation with modified recursive convolution scheme on basis of
            ! second-order accurate convolution term calculation from equation (21) of
            ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
            ! Anisotropic-medium PML for vector FETD with modified basis functions,
            ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)
            if (time_stepping_scheme == 1) then
              ! Newmark
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

            ! update e1, e11, e13 in ADE formation with fourth-order LDDRK scheme
            else if (time_stepping_scheme == 2) then
              ! LDDRK
              e1_LDDRK(i,j,ispec,i_sls) = ALPHA_LDDRK(i_stage) * e1_LDDRK(i,j,ispec,i_sls) + &
                                          deltat * (theta_n_u * phinu1 - e1(i,j,ispec,i_sls) * tauinvnu1)
              e1(i,j,ispec,i_sls) = e1(i,j,ispec,i_sls) + BETA_LDDRK(i_stage) * e1_LDDRK(i,j,ispec,i_sls)

              e11_LDDRK(i,j,ispec,i_sls) = ALPHA_LDDRK(i_stage) * e11_LDDRK(i,j,ispec,i_sls) + &
                                           deltat * ((dux_dxl_n(i,j,ispec)-theta_n_u/TWO) * phinu2) - &
                                           deltat * (e11(i,j,ispec,i_sls) * tauinvnu2)
              e11(i,j,ispec,i_sls) = e11(i,j,ispec,i_sls)+BETA_LDDRK(i_stage)*e11_LDDRK(i,j,ispec,i_sls)

              e13_LDDRK(i,j,ispec,i_sls) = ALPHA_LDDRK(i_stage) * e13_LDDRK(i,j,ispec,i_sls) + &
                                           deltat * ((dux_dzl_n(i,j,ispec) + duz_dxl_n(i,j,ispec))*phinu2) - &
                                           deltat * (e13(i,j,ispec,i_sls) * tauinvnu2)
              e13(i,j,ispec,i_sls) = e13(i,j,ispec,i_sls)+BETA_LDDRK(i_stage) * e13_LDDRK(i,j,ispec,i_sls)

            ! update e1, e11, e13 in ADE formation with classical fourth-order Runge-Kutta scheme
            else if (time_stepping_scheme == 3) then
              ! RK
              e1_force_RK(i,j,ispec,i_sls,i_stage) = deltat * (theta_n_u * phinu1 - e1(i,j,ispec,i_sls) * tauinvnu1)

              if (i_stage == 1 .or. i_stage == 2 .or. i_stage == 3) then
                if (i_stage == 1)weight_rk = 0.5_CUSTOM_REAL
                if (i_stage == 2)weight_rk = 0.5_CUSTOM_REAL
                if (i_stage == 3)weight_rk = 1._CUSTOM_REAL

                if (i_stage == 1) e1_initial_rk(i,j,ispec,i_sls) = e1(i,j,ispec,i_sls)
                e1(i,j,ispec,i_sls) = e1_initial_rk(i,j,ispec,i_sls) + weight_rk * e1_force_RK(i,j,ispec,i_sls,i_stage)
              else if (i_stage == 4) then
                e1(i,j,ispec,i_sls) = e1_initial_rk(i,j,ispec,i_sls) + 1._CUSTOM_REAL / 6._CUSTOM_REAL * &
                                      (e1_force_RK(i,j,ispec,i_sls,1) + 2._CUSTOM_REAL * e1_force_RK(i,j,ispec,i_sls,2) + &
                                       2._CUSTOM_REAL * e1_force_RK(i,j,ispec,i_sls,3) + e1_force_RK(i,j,ispec,i_sls,4))
              endif

              e11_force_RK(i,j,ispec,i_sls,i_stage) = deltat * ((dux_dxl_n(i,j,ispec)-theta_n_u/TWO) * phinu2 - &
                                                                 e11(i,j,ispec,i_sls) * tauinvnu2)

              if (i_stage == 1 .or. i_stage == 2 .or. i_stage == 3) then
                if (i_stage == 1)weight_rk = 0.5_CUSTOM_REAL
                if (i_stage == 2)weight_rk = 0.5_CUSTOM_REAL
                if (i_stage == 3)weight_rk = 1._CUSTOM_REAL

                if (i_stage == 1) e11_initial_rk(i,j,ispec,i_sls) = e11(i,j,ispec,i_sls)
                e11(i,j,ispec,i_sls) = e11_initial_rk(i,j,ispec,i_sls) + weight_rk * e11_force_RK(i,j,ispec,i_sls,i_stage)
              else if (i_stage == 4) then
                e11(i,j,ispec,i_sls) = e11_initial_rk(i,j,ispec,i_sls) + 1._CUSTOM_REAL / 6._CUSTOM_REAL * &
                                       (e11_force_RK(i,j,ispec,i_sls,1) + 2._CUSTOM_REAL * e11_force_RK(i,j,ispec,i_sls,2) + &
                                        2._CUSTOM_REAL * e11_force_RK(i,j,ispec,i_sls,3) + e11_force_RK(i,j,ispec,i_sls,4))
              endif

              e13_force_RK(i,j,ispec,i_sls,i_stage) = deltat * ((dux_dzl_n(i,j,ispec) + duz_dxl_n(i,j,ispec))*phinu2 - &
                                                                 e13(i,j,ispec,i_sls) * tauinvnu2)
              if (i_stage == 1 .or. i_stage == 2 .or. i_stage == 3) then
                if (i_stage == 1)weight_rk = 0.5_CUSTOM_REAL
                if (i_stage == 2)weight_rk = 0.5_CUSTOM_REAL
                if (i_stage == 3)weight_rk = 1._CUSTOM_REAL

                if (i_stage == 1) e13_initial_rk(i,j,ispec,i_sls) = e13(i,j,ispec,i_sls)
                e13(i,j,ispec,i_sls) = e13_initial_rk(i,j,ispec,i_sls) + weight_rk * e13_force_RK(i,j,ispec,i_sls,i_stage)
              else if (i_stage == 4) then
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
  endif ! iphase
!!!! end of update memory variable in viscoelastic simulation

! this to avoid a warning at execution time about an undefined variable being used
! for the SH component in the case of a P-SV calculation, and vice versa
  sigma_xx = 0._CUSTOM_REAL
  sigma_xz = 0._CUSTOM_REAL
  sigma_zz = 0._CUSTOM_REAL
  sigma_zx = 0._CUSTOM_REAL
  ! SH-case
  sigma_xy = 0._CUSTOM_REAL
  sigma_zy = 0._CUSTOM_REAL

  ! choses inner/outer elements
  if (iphase == 1) then
    num_elements = nspec_outer_elastic
  else
    num_elements = nspec_inner_elastic
  endif

  ! loop over spectral elements
  do ispec_p = 1,num_elements

    ! returns element id from stored element list
    ispec = phase_ispec_inner_elastic(ispec_p,iphase)

    tempx1(:,:) = 0._CUSTOM_REAL
    tempz1(:,:) = 0._CUSTOM_REAL
    tempx2(:,:) = 0._CUSTOM_REAL
    tempz2(:,:) = 0._CUSTOM_REAL
    ! AXISYM-case
    tempx3(:,:) = 0._CUSTOM_REAL

    sigma_thetatheta(:,:) = 0._CUSTOM_REAL

    !--- elastic spectral element
    if (ispec_is_elastic(ispec)) then
      ! get unrelaxed elastic parameters of current spectral element
      lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
      mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
      lambdalplusmul_unrelaxed_elastic = lambdal_unrelaxed_elastic + mul_unrelaxed_elastic
      lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec))

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX
          !--- if external medium, get elastic parameters of current grid point
          if (assign_external_model) then
            cpl = vpext(i,j,ispec)
            csl = vsext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
            mul_unrelaxed_elastic = rhol*csl*csl
            lambdal_unrelaxed_elastic = rhol*cpl*cpl - TWO*mul_unrelaxed_elastic
            lambdalplusmul_unrelaxed_elastic = lambdal_unrelaxed_elastic + mul_unrelaxed_elastic
            lambdaplus2mu_unrelaxed_elastic = lambdal_unrelaxed_elastic + TWO*mul_unrelaxed_elastic
          endif

          ! derivative along x and along z
          dux_dxi = 0._CUSTOM_REAL
          duz_dxi = 0._CUSTOM_REAL
          dux_dgamma = 0._CUSTOM_REAL
          duz_dgamma = 0._CUSTOM_REAL

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
            if (AXISYM) then
              if (is_on_the_axis(ispec)) then
                do k = 1,NGLJ
                  dux_dxi = dux_dxi + b_displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                  duz_dxi = duz_dxi + b_displ_elastic(2,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                  dux_dgamma = dux_dgamma + b_displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
                  duz_dgamma = duz_dgamma + b_displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
                enddo
              else
                do k = 1,NGLJ
                  dux_dxi = dux_dxi + b_displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
                  duz_dxi = duz_dxi + b_displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
                  dux_dgamma = dux_dgamma + b_displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
                  duz_dgamma = duz_dgamma + b_displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
                enddo
              endif
            else
              do k = 1,NGLLX
                dux_dxi = dux_dxi + b_displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
                duz_dxi = duz_dxi + b_displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
                dux_dgamma = dux_dgamma + b_displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
                duz_dgamma = duz_dgamma + b_displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
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

          if (AXISYM .and. is_on_the_axis(ispec) .and. i == 1) then ! d_uz/dr=0 on the axis
            duz_dxl = 0.d0
          endif

          ! compute stress tensor (include attenuation or anisotropy if needed)
          if (ATTENUATION_VISCOELASTIC) then
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

            if (AXISYM) then
              if (is_on_the_axis(ispec)) then
                if (is_on_the_axis(ispec) .and. i == 1) then ! First GLJ point
                  sigma_xx = 0._CUSTOM_REAL
                  sigma_zz = 0._CUSTOM_REAL
                  sigma_thetatheta(i,j) = 0._CUSTOM_REAL
                  xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                  r_xiplus1(i,j) = xxi
                  do k = 1,NGLJ
                    sigma_xx = sigma_xx + b_displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                    sigma_zz = sigma_zz + b_displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                    sigma_thetatheta(i,j) = sigma_thetatheta(i,j) + b_displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
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
                             + lambdal_unrelaxed_elastic*b_displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                  sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                             + lambdal_unrelaxed_elastic*b_displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                  sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
                  sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                          + lambdaplus2mu_unrelaxed_elastic &
                                          * b_displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                  r_xiplus1(i,j) = coord(1,ibool(i,j,ispec))/(xiglj(i)+ONE)
                endif
              else ! Not on the axis
                sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl &
                           + lambdal_unrelaxed_elastic*b_displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                           + lambdal_unrelaxed_elastic*b_displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
                sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                        + lambdaplus2mu_unrelaxed_elastic &
                                        * b_displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
              endif
            else ! Not axisym
              sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
              sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
              sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl
            endif

            ! add the memory variables (Carcione 2007 page 125)
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

            sigma_xx = sigma_xx + lambdalplusmul_unrelaxed_elastic * e1_sum + TWO * mul_unrelaxed_elastic * e11_sum
            sigma_xz = sigma_xz + mul_unrelaxed_elastic * e13_sum
            sigma_zz = sigma_zz + lambdalplusmul_unrelaxed_elastic * e1_sum - TWO * mul_unrelaxed_elastic * e11_sum
            sigma_zx = sigma_xz

          else
            ! no attenuation

            if (AXISYM) then
              if (is_on_the_axis(ispec)) then
                if (is_on_the_axis(ispec) .and. i == 1) then ! First GLJ point
                  sigma_xx = 0._CUSTOM_REAL
                  sigma_zz = 0._CUSTOM_REAL
                  sigma_thetatheta(i,j) = 0._CUSTOM_REAL
                  xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                  r_xiplus1(i,j) = xxi
                  do k = 1,NGLJ
                    sigma_xx = sigma_xx + b_displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                    sigma_zz = sigma_zz + b_displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                    sigma_thetatheta(i,j) = sigma_thetatheta(i,j) + b_displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
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
                             + lambdal_unrelaxed_elastic*b_displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                  sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                             + lambdal_unrelaxed_elastic*b_displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                  sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
                  sigma_zx = sigma_xz
                  sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                          + lambdaplus2mu_unrelaxed_elastic &
                                          * b_displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                  r_xiplus1(i,j) = coord(1,ibool(i,j,ispec))/(xiglj(i)+ONE)
                endif
              else ! Not on the axis
                sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl &
                           + lambdal_unrelaxed_elastic*b_displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                           + lambdal_unrelaxed_elastic*b_displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
                sigma_zx = sigma_xz
                sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                        + lambdaplus2mu_unrelaxed_elastic &
                                        * b_displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
              endif
            else
              ! Not axisym
              if (P_SV) then
                ! P_SV-case
                sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
                sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
                sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl
                sigma_zx = sigma_xz
              else
                ! SH-case
                sigma_xy = mul_unrelaxed_elastic*dux_dxl
                sigma_zy = mul_unrelaxed_elastic*dux_dzl
              endif
            endif
          endif

          ! full anisotropy
          if (ispec_is_anisotropic(ispec)) then
            if (assign_external_model) then
              c11 = c11ext(i,j,ispec)
              c13 = c13ext(i,j,ispec)
              c15 = c15ext(i,j,ispec)
              c33 = c33ext(i,j,ispec)
              c35 = c35ext(i,j,ispec)
              c55 = c55ext(i,j,ispec)
              c12 = c12ext(i,j,ispec)
              c23 = c23ext(i,j,ispec)
              c25 = c25ext(i,j,ispec)
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
            endif

            ! implement anisotropy in 2D
            sigma_xx = c11*dux_dxl + c13*duz_dzl + c15*(duz_dxl + dux_dzl)
            sigma_zz = c13*dux_dxl + c33*duz_dzl + c35*(duz_dxl + dux_dzl)
            sigma_xz = c15*dux_dxl + c35*duz_dzl + c55*(duz_dxl + dux_dzl)
            sigma_zx = sigma_xz  !ZN I add this line, since no where compute the sigma_zx for anistropic simulation
          endif

          ! weak formulation term based on stress tensor (non-symmetric form)
          ! also add GLL integration weights
          jacobianl = jacobian(i,j,ispec)

          !! ABAB with the notations of Komatitsch & Tromp 1999 (with 3 -> 2) :
          ! tempx1(i,j) = w.J.F_{11}^{ij}
          ! tempz1(i,j) = w.J.F_{21}^{ij}
          ! tempx2(i,j) = w.J.F_{12}^{ij}
          ! tempz2(i,j) = w.J.F_{22}^{ij}

          if (AXISYM) then
            if (is_on_the_axis(ispec)) then
              tempx3(i,j) = wzgll(j)*jacobian(1,j,ispec)*sigma_thetatheta(1,j)*hprimeBarwglj_xx(1,i)

              if (abs(coord(1,ibool(i,j,ispec))) > TINYVAL) then ! Not first GLJ point
                if (i == 1) then
                  call stop_the_code("error: an axial element is rotated. The code should have been stopped before. &
                   &Check that your coordinates are greater than TINYVAL. Maybe you should also have a look to &
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
            if (P_SV) then
              ! P_SV-case
              tempx1(i,j) = wzgll(j)*jacobianl*(sigma_xx*xixl+sigma_zx*xizl) ! this goes to accel_x
              tempz1(i,j) = wzgll(j)*jacobianl*(sigma_xz*xixl+sigma_zz*xizl) ! this goes to accel_z

              tempx2(i,j) = wxgll(i)*jacobianl*(sigma_xx*gammaxl+sigma_zx*gammazl) ! this goes to accel_x
              tempz2(i,j) = wxgll(i)*jacobianl*(sigma_xz*gammaxl+sigma_zz*gammazl) ! this goes to accel_z
            else
              ! SH-case
              tempx1(i,j) = wzgll(j)*jacobianl*(sigma_xy*xixl+sigma_zy*xizl) ! this goes to accel_x
              tempx2(i,j) = wxgll(i)*jacobianl*(sigma_xy*gammaxl+sigma_zy*gammazl) ! this goes to accel_x
              tempz1(i,j) = 0._CUSTOM_REAL
              tempz2(i,j) = 0._CUSTOM_REAL
            endif
          endif
        enddo
      enddo  ! end of the loops on the collocation points i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! second double-loop over GLL to compute all the terms
      !
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          ! along x direction and z direction
          ! and assemble the contributions
          ! we can merge the two loops because NGLLX == NGLLZ
          if (AXISYM) then
            if (is_on_the_axis(ispec)) then
              do k = 1,NGLJ
                b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) &
                                         - (tempx1(k,j)*hprimeBarwglj_xx(k,i) + tempx2(i,k)*hprimewgll_zz(k,j))
                b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) &
                                         - (tempz1(k,j)*hprimeBarwglj_xx(k,i) + tempz2(i,k)*hprimewgll_zz(k,j))
              enddo
              b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - tempx3(i,j)
            else ! Axisym but not on the axis
              do k = 1,NGLLX
                b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) &
                                         - (tempx1(k,j)*hprimewgll_xx(k,i) + tempx2(i,k)*hprimewgll_zz(k,j))
                b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) &
                                         - (tempz1(k,j)*hprimewgll_xx(k,i) + tempz2(i,k)*hprimewgll_zz(k,j))
              enddo
              b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - tempx3(i,j)
            endif
          else !if AXISYM == false
            do k = 1,NGLLX
              b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - &
                                         (tempx1(k,j)*hprimewgll_xx(k,i) + tempx2(i,k)*hprimewgll_zz(k,j))
              b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) - &
                                         (tempz1(k,j)*hprimewgll_xx(k,i) + tempz2(i,k)*hprimewgll_zz(k,j))
            enddo
          endif

        enddo
      enddo ! second loop over the GLL points
    endif ! end of test if elastic element
  enddo ! end of loop over all spectral elements

  end subroutine compute_forces_viscoelastic_backward

