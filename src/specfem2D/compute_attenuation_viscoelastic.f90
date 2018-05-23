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

! for viscoelastic solver

  subroutine compute_attenuation_viscoelastic(e1,e11,e13,dux_dxl,dux_dzl,duz_dxl,duz_dzl,dux_dxl_old,duz_dzl_old, &
                                              dux_dzl_plus_duz_dxl_old,PML_BOUNDARY_CONDITIONS,i,j,ispec, &
                                              e1_sum,e11_sum,e13_sum)

  ! updates memory variable in viscoelastic simulation

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,TWO,ALPHA_LDDRK,BETA_LDDRK,C_LDDRK

  use specfem_par, only: nspec_ATT_el,ATTENUATION_VISCOELASTIC,N_SLS, &
                         ispec_is_PML, &
                         inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2,time_stepping_scheme,i_stage,deltat, &
                         e1_LDDRK,e11_LDDRK,e13_LDDRK,e1_initial_rk,e11_initial_rk,e13_initial_rk, &
                         e1_force_RK,e11_force_RK,e13_force_RK,A_newmark_nu1,B_newmark_nu1,A_newmark_nu2,B_newmark_nu2

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: dux_dxl,dux_dzl,duz_dxl,duz_dzl
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ,nspec_ATT_el),intent(inout) :: dux_dxl_old,duz_dzl_old,dux_dzl_plus_duz_dxl_old
  real(kind=CUSTOM_REAL),intent(out) :: e1_sum,e11_sum,e13_sum
  integer, intent(in) :: i,j,ispec
  ! CPML coefficients and memory variables
  logical, intent(in) :: PML_BOUNDARY_CONDITIONS
  real(kind=CUSTOM_REAL),dimension(N_SLS,NGLLX,NGLLZ,nspec_ATT_el),intent(inout) :: e1,e11,e13

  ! local variables
  integer :: i_sls

  ! for attenuation
  real(kind=CUSTOM_REAL) :: a_newmark,phinu1,phinu2,theta_n_u,theta_nsub1_u
  double precision :: tauinvnu1,tauinvnu2

  ! temporary RK4 variable
  real(kind=CUSTOM_REAL) :: weight_rk

! update the memory variables using a convolution or using a differential equation
! (tests made by Ting Yu and also by Zhinan Xie, CNRS Marseille, France, show that it is better to leave it to .true.)
  logical, parameter :: CONVOLUTION_MEMORY_VARIABLES = .true.

  e1_sum = 0.0
  e11_sum = 0.0
  e13_sum = 0.0

  ! checks if anything to do
  if (.not. ATTENUATION_VISCOELASTIC) return
  if (PML_BOUNDARY_CONDITIONS .and. ispec_is_PML(ispec)) return
  ! convention to indicate that Q = 9999 in that element i.e. that there is no viscoelasticity at that GLL point
  if (inv_tau_sigma_nu1(i,j,ispec,1) < 0.) return

  if (.not. CONVOLUTION_MEMORY_VARIABLES) &
    call stop_the_code('CONVOLUTION_MEMORY_VARIABLES == .false. is not accurate enough and has been discontinued for now')

  theta_n_u = dux_dxl + duz_dzl
  if (time_stepping_scheme == 1) theta_nsub1_u = dux_dxl_old(i,j,ispec) + duz_dzl_old(i,j,ispec)

  ! loop on all the standard linear solids
  do i_sls = 1,N_SLS

    if (time_stepping_scheme /= 1) then
      phinu1 = phi_nu1(i,j,ispec,i_sls)
      tauinvnu1 = inv_tau_sigma_nu1(i,j,ispec,i_sls)
      phinu2 = phi_nu2(i,j,ispec,i_sls)
      tauinvnu2 = inv_tau_sigma_nu2(i,j,ispec,i_sls)
    endif

    ! update e1, e11, e13 in convolution formulation with modified recursive convolution scheme on basis of
    ! second-order accurate convolution term calculation from equation (21) of
    ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
    ! Anisotropic-medium PML for vector FETD with modified basis functions,
    ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)
    select case (time_stepping_scheme)
    case (1)
      ! Newmark

! update the memory variables using a convolution or using a differential equation
! From Zhinan Xie and Dimitri Komatitsch:
! For cases in which a value of tau_sigma is small, then its inverse is large,
! which may result in a in stiff ordinary differential equation to solve;
! in such a case, resorting to the convolution formulation is better.
!         if (CONVOLUTION_MEMORY_VARIABLES) then
!! DK DK inlined this for speed            call compute_coef_convolution(tauinvnu1,deltat,coef0,coef1,coef2)
        a_newmark = A_newmark_nu1(i_sls,i,j,ispec)
        e1(i_sls,i,j,ispec) = a_newmark * a_newmark * e1(i_sls,i,j,ispec) + &
                              B_newmark_nu1(i_sls,i,j,ispec) * (theta_n_u + a_newmark * theta_nsub1_u)

!! DK DK inlined this for speed            call compute_coef_convolution(tauinvnu2,deltat,coef0,coef1,coef2)
        a_newmark = A_newmark_nu2(i_sls,i,j,ispec)
        e11(i_sls,i,j,ispec) = a_newmark * a_newmark * e11(i_sls,i,j,ispec) + &
                               B_newmark_nu2(i_sls,i,j,ispec) * (dux_dxl-theta_n_u/TWO + &
                               a_newmark * (dux_dxl_old(i,j,ispec)-theta_nsub1_u/TWO))

        e13(i_sls,i,j,ispec) = a_newmark * a_newmark * e13(i_sls,i,j,ispec) + &
                               B_newmark_nu2(i_sls,i,j,ispec) * (dux_dzl + duz_dxl + &
                               a_newmark * (dux_dzl_plus_duz_dxl_old(i,j,ispec)))

!         else
!           stop 'CONVOLUTION_MEMORY_VARIABLES == .false. is not accurate enough and has been discontinued for now'
!           e1(i_sls,i,j,ispec) = e1(i_sls,i,j,ispec) + deltat * &
!                (- e1(i_sls,i,j,ispec)*tauinvnu1 + phinu1 * theta_n_u)
!
!           e11(i_sls,i,j,ispec) = e11(i_sls,i,j,ispec) + deltat * &
!                (- e11(i_sls,i,j,ispec)*tauinvnu2 + phinu2 * (dux_dxl_n(i,j,ispec)-theta_n_u/TWO))
!
!           e13(i_sls,i,j,ispec) = e13(i_sls,i,j,ispec) + deltat * &
!              (- e13(i_sls,i,j,ispec)*tauinvnu2 + phinu2 * (dux_dzl_n(i,j,ispec) + duz_dxl_n(i,j,ispec)))
!         endif

    case (2)
      ! LDDRK
      ! update e1, e11, e13 in ADE formation with fourth-order LDDRK scheme
      e1_LDDRK(i,j,ispec,i_sls) = ALPHA_LDDRK(i_stage) * e1_LDDRK(i,j,ispec,i_sls) + &
                                  deltat * (theta_n_u * phinu1 - e1(i_sls,i,j,ispec) * tauinvnu1)
      e1(i_sls,i,j,ispec) = e1(i_sls,i,j,ispec) + BETA_LDDRK(i_stage) * e1_LDDRK(i,j,ispec,i_sls)

      e11_LDDRK(i,j,ispec,i_sls) = ALPHA_LDDRK(i_stage) * e11_LDDRK(i,j,ispec,i_sls) + &
                                   deltat * ((dux_dxl-theta_n_u/TWO) * phinu2) - &
                                   deltat * (e11(i_sls,i,j,ispec) * tauinvnu2)
      e11(i_sls,i,j,ispec) = e11(i_sls,i,j,ispec)+BETA_LDDRK(i_stage)*e11_LDDRK(i,j,ispec,i_sls)

      e13_LDDRK(i,j,ispec,i_sls) = ALPHA_LDDRK(i_stage) * e13_LDDRK(i,j,ispec,i_sls) + &
                                   deltat * ((dux_dzl + duz_dxl)*phinu2) - &
                                   deltat * (e13(i_sls,i,j,ispec) * tauinvnu2)
      e13(i_sls,i,j,ispec) = e13(i_sls,i,j,ispec)+BETA_LDDRK(i_stage) * e13_LDDRK(i,j,ispec,i_sls)

    case (3)
      ! Runge-Kutta
      ! update e1, e11, e13 in ADE formation with classical fourth-order Runge-Kutta scheme
      e1_force_RK(i,j,ispec,i_sls,i_stage) = deltat * (theta_n_u * phinu1 - e1(i_sls,i,j,ispec) * tauinvnu1)

      if (i_stage == 1 .or. i_stage == 2 .or. i_stage == 3) then
        if (i_stage == 1) weight_rk = 0.5_CUSTOM_REAL
        if (i_stage == 2) weight_rk = 0.5_CUSTOM_REAL
        if (i_stage == 3) weight_rk = 1._CUSTOM_REAL
        if (i_stage == 1) e1_initial_rk(i,j,ispec,i_sls) = e1(i_sls,i,j,ispec)
        e1(i_sls,i,j,ispec) = e1_initial_rk(i,j,ispec,i_sls) + weight_rk * e1_force_RK(i,j,ispec,i_sls,i_stage)
      else if (i_stage == 4) then
        e1(i_sls,i,j,ispec) = e1_initial_rk(i,j,ispec,i_sls) + 1._CUSTOM_REAL / 6._CUSTOM_REAL * &
                              (e1_force_RK(i,j,ispec,i_sls,1) + 2._CUSTOM_REAL * e1_force_RK(i,j,ispec,i_sls,2) + &
                               2._CUSTOM_REAL * e1_force_RK(i,j,ispec,i_sls,3) + e1_force_RK(i,j,ispec,i_sls,4))
      endif

      e11_force_RK(i,j,ispec,i_sls,i_stage) = deltat * ((dux_dxl-theta_n_u/TWO) * phinu2 - &
                                                         e11(i_sls,i,j,ispec) * tauinvnu2)

      if (i_stage == 1 .or. i_stage == 2 .or. i_stage == 3) then
        if (i_stage == 1) weight_rk = 0.5_CUSTOM_REAL
        if (i_stage == 2) weight_rk = 0.5_CUSTOM_REAL
        if (i_stage == 3) weight_rk = 1._CUSTOM_REAL

        if (i_stage == 1) e11_initial_rk(i,j,ispec,i_sls) = e11(i_sls,i,j,ispec)
        e11(i_sls,i,j,ispec) = e11_initial_rk(i,j,ispec,i_sls) + weight_rk * e11_force_RK(i,j,ispec,i_sls,i_stage)
      else if (i_stage == 4) then
        e11(i_sls,i,j,ispec) = e11_initial_rk(i,j,ispec,i_sls) + 1._CUSTOM_REAL / 6._CUSTOM_REAL * &
                               (e11_force_RK(i,j,ispec,i_sls,1) + 2._CUSTOM_REAL * e11_force_RK(i,j,ispec,i_sls,2) + &
                                2._CUSTOM_REAL * e11_force_RK(i,j,ispec,i_sls,3) + e11_force_RK(i,j,ispec,i_sls,4))
      endif

      e13_force_RK(i,j,ispec,i_sls,i_stage) = deltat * ((dux_dzl + duz_dxl)*phinu2 - &
                                                         e13(i_sls,i,j,ispec) * tauinvnu2)
      if (i_stage == 1 .or. i_stage == 2 .or. i_stage == 3) then
        if (i_stage == 1) weight_rk = 0.5_CUSTOM_REAL
        if (i_stage == 2) weight_rk = 0.5_CUSTOM_REAL
        if (i_stage == 3) weight_rk = 1._CUSTOM_REAL

        if (i_stage == 1) e13_initial_rk(i,j,ispec,i_sls) = e13(i_sls,i,j,ispec)
        e13(i_sls,i,j,ispec) = e13_initial_rk(i,j,ispec,i_sls) + weight_rk * e13_force_RK(i,j,ispec,i_sls,i_stage)
      else if (i_stage == 4) then
        e13(i_sls,i,j,ispec) = e13_initial_rk(i,j,ispec,i_sls) + 1._CUSTOM_REAL / 6._CUSTOM_REAL * &
                               (e13_force_RK(i,j,ispec,i_sls,1) + 2._CUSTOM_REAL * e13_force_RK(i,j,ispec,i_sls,2) + &
                                2._CUSTOM_REAL * e13_force_RK(i,j,ispec,i_sls,3) + e13_force_RK(i,j,ispec,i_sls,4))
      endif

    case default
      call stop_the_code('Time stepping scheme not implemented yet in viscoelastic attenuation update')
    end select

  enddo ! i_sls

  e1_sum = sum(e1(:,i,j,ispec))
  e11_sum = sum(e11(:,i,j,ispec))
  e13_sum = sum(e13(:,i,j,ispec))

  if (time_stepping_scheme == 1) then
        !Update of grad(Displ)
        dux_dxl_old(i,j,ispec) = dux_dxl
        duz_dzl_old(i,j,ispec) = duz_dzl
        dux_dzl_plus_duz_dxl_old(i,j,ispec) = dux_dzl + duz_dxl
  endif

  end subroutine compute_attenuation_viscoelastic

