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

  subroutine compute_coef_convolution(bb,DT,coef0,coef1,coef2)

  ! compute coefficient used in second-order convolution scheme, from
  ! second-order accurate convolution term calculation from equation (21) of
  ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
  ! Anisotropic-medium PML for vector FETD with modified basis functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  implicit none

  double precision :: bb,coef0,coef1,coef2,temp
  double precision :: DT

  ! recursive convolution coefficients
  !
  ! see Xie et al. (2014), second-order recursive scheme given by eq. (60)
  !                        and also appendix D, eq. (D6a) and (D6b) for p = 0
  !
  ! coefficients needed for the recursive scheme are:
  !    coef0 = exp(-b delta_t)
  !          = exp(- 1/2 b delta_t) * exp(- 1/2 b delta_t)
  !
  !    coef1 = 1/b (1 - exp( - 1/2 b delta_t)                             -> see also factor xi_0^(n+1) in eq. D6b
  !
  !    coef2 = 1/b (1 - exp( - 1/2 b delta_t) exp(- 1/2 b delta_t)
  !          = coef1 * exp(- 1/2 b delta_t)                               -> see also factor xi_0^n in eq. D6a
  !
  ! helper variables
  temp = exp(- 0.5d0 * bb * DT)

  coef0 = temp*temp

  if (abs(bb) > 1d-5) then
    ! second-order convolution scheme coefficients
    coef1 = (1.d0 - temp) / bb
    coef2 = coef1 * temp

    ! first-order scheme
    !coef1 = (1.d0 - coef0) / bb
    !coef2 = 0.d0
  else
    ! approximation for small beta
    coef1 = 0.5d0 * DT
    coef2 = coef1

    ! Taylor expansion to third-order
    !coef1 = 0.5d0 * DT + &
    !        (- 1.d0/8.d0 * deltatpow2 * bb + 1.d0/48.d0 * DT**3 * bb**2 - 1.d0/384.d0 * DT**4 * bb**3)
    !coef2 = 0.5d0 * DT + &
    !        (- 3.d0/8.d0 * deltatpow2 * bb + 7.d0/48.d0 * DT**3 * bb**2 - 5.d0/128.d0 * DT**4 * bb**3)
  endif

  end subroutine compute_coef_convolution

!========================================================================

  subroutine lik_parameter_computation(DT,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                       CPML_region_local,index_ik,A_0,A_1,A_2,bb_1,bb_2, &
                                       coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2)

  use constants, only: CPML_X_ONLY,CPML_Z_ONLY,CPML_XZ,CUSTOM_REAL

  implicit none

  double precision, intent(in) :: DT
  double precision, intent(in) :: kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z
  integer, intent(in) :: CPML_region_local,index_ik

  double precision, intent(out) :: A_0,A_1,A_2
  double precision, intent(out) :: coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2

  ! local variables
  double precision :: bar_A_0,bar_A_1,bar_A_2,bb_1,bb_2,gamma_x,gamma_z
  integer :: CPML_X_ONLY_TEMP,CPML_Z_ONLY_TEMP,CPML_XZ_TEMP

  logical,parameter :: FIRST_ORDER_CONVOLUTION = .false.

  if (index_ik == 13) then
    CPML_X_ONLY_TEMP = CPML_X_ONLY
    CPML_Z_ONLY_TEMP = CPML_Z_ONLY
    CPML_XZ_TEMP = CPML_XZ
  else if (index_ik == 31) then
    CPML_X_ONLY_TEMP = CPML_Z_ONLY
    CPML_Z_ONLY_TEMP = CPML_X_ONLY
    CPML_XZ_TEMP = CPML_XZ
  else
    call stop_the_code('In lik_parameter_computation index_ik must be equal to 13 or 31')
  endif

  if (CPML_region_local == CPML_XZ_TEMP) then
    !----------------A0-------------------------
    bar_A_0 = kappa_x / kappa_z
    A_0 = bar_A_0
    gamma_x = (alpha_x * beta_z + alpha_x**2 + 2._CUSTOM_REAL * beta_x * alpha_z - &
               2._CUSTOM_REAL * alpha_x * (beta_x + alpha_z)) / (beta_z - alpha_x)
    gamma_z = (alpha_x * beta_z + beta_z**2 + 2._CUSTOM_REAL * beta_x * alpha_z - &
               2._CUSTOM_REAL * beta_z * (beta_x + alpha_z)) / (alpha_x - beta_z)

    !----------------A1,2-------------------------
    bar_A_1 = bar_A_0 / 2._CUSTOM_REAL * (gamma_x - alpha_x)
    bar_A_2 = bar_A_0 / 2._CUSTOM_REAL * (gamma_z - beta_z)

    A_1 = bar_A_1
    A_2 = bar_A_2
  else if (CPML_region_local == CPML_X_ONLY_TEMP) then
    !----------------A0-------------------------
    bar_A_0 = kappa_x
    A_0 = bar_A_0
    !----------------A1,2,3-------------------------
    bar_A_1 = - bar_A_0 * (alpha_x - beta_x)
    bar_A_2 = 0.d0

    A_1 = bar_A_1
    A_2 = bar_A_2
  else if (CPML_region_local == CPML_Z_ONLY_TEMP) then
    !----------------A0-------------------------
    bar_A_0 = 1.d0 / kappa_z
    A_0 = bar_A_0

    !----------------A1,2,3-------------------------
    bar_A_1 = 0.d0
    bar_A_2 = - bar_A_0 * (beta_z - alpha_z)

    A_1 = bar_A_1
    A_2 = bar_A_2
  endif

  ! gets recursive convolution coefficients
  ! alpha coefficients
  bb_1 = alpha_x
  call compute_coef_convolution(bb_1,DT,coef0_1,coef1_1,coef2_1)
  ! beta coefficients
  bb_2 = beta_z
  call compute_coef_convolution(bb_2,DT,coef0_2,coef1_2,coef2_2)

 end subroutine lik_parameter_computation

!========================================================================

  subroutine lik_parameter_computation_viscoelastic(kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                                    CPML_region_local,index_ik,A_0,A_1,A_2,A_zener, &
                                                    inv_tau_temp,tau_epsilon_temp)

  use constants, only: NDIM,CUSTOM_REAL,CPML_X_ONLY,CPML_Z_ONLY,CPML_XZ

  implicit none

  double precision, intent(in) :: kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z
  integer, intent(in) :: CPML_region_local,index_ik
  double precision, intent(in) :: tau_epsilon_temp,inv_tau_temp

  double precision, intent(out) :: A_0,A_1,A_2,A_zener

  ! local variables
  integer :: i,CPML_X_ONLY_TEMP,CPML_Z_ONLY_TEMP,CPML_XZ_TEMP,N_PMLSF
  double precision :: bar_A_0,bar_A_1,bar_A_2
  double precision, dimension(NDIM) :: beta_xz,alpha_xz,gamma_xz,eta_xz
  double precision :: gamma_only,alpha_only,beta_only,eta_only


  if (index_ik == 13) then
    CPML_X_ONLY_TEMP = CPML_X_ONLY
    CPML_Z_ONLY_TEMP = CPML_Z_ONLY
    CPML_XZ_TEMP = CPML_XZ
  else if (index_ik == 31) then
    CPML_X_ONLY_TEMP = CPML_Z_ONLY
    CPML_Z_ONLY_TEMP = CPML_X_ONLY
    CPML_XZ_TEMP = CPML_XZ
  else
    call stop_the_code('In lik_parameter_computation index_ik must be equal to 13 or 31')
  endif

  if (CPML_region_local == CPML_XZ_TEMP) then

    bar_A_0 = kappa_x / kappa_z
    N_PMLSF = 2
    alpha_xz(1) = alpha_x
    alpha_xz(2) = beta_z
    beta_xz(1) = beta_x
    beta_xz(2) = alpha_z

    call decompose_rational_fraction_PML(N_PMLSF,beta_xz,alpha_xz)

    do i = 1,N_PMLSF
      gamma_xz(i) = (beta_xz(i) * (inv_tau_temp - alpha_xz(i) * tau_epsilon_temp * inv_tau_temp) + &
                     alpha_xz(i) * alpha_xz(i) * (tau_epsilon_temp * inv_tau_temp - 1._CUSTOM_REAL)) / &
                    (inv_tau_temp - alpha_xz(i))
      eta_xz(i) = inv_tau_temp * (beta_xz(i) - alpha_xz(i)) / (inv_tau_temp - alpha_xz(i))
    enddo

    bar_A_1 = bar_A_0 / 2._CUSTOM_REAL * (gamma_xz(1) - alpha_xz(1))
    bar_A_2 = bar_A_0 / 2._CUSTOM_REAL * (gamma_xz(2) - alpha_xz(2))

    A_0 = bar_A_0
    A_1 = bar_A_1
    A_2 = bar_A_2

    A_zener = bar_A_0 * (tau_epsilon_temp * inv_tau_temp - 1._CUSTOM_REAL) * &
              ((eta_xz(1) + eta_xz(2)) / 2._CUSTOM_REAL  - inv_tau_temp)

  else if (CPML_region_local == CPML_X_ONLY_TEMP) then

    bar_A_0 = kappa_x
    alpha_only = alpha_x
    beta_only = beta_x

    gamma_only = (beta_only * (inv_tau_temp - alpha_only * tau_epsilon_temp * inv_tau_temp) + &
                  alpha_only * alpha_only * (tau_epsilon_temp * inv_tau_temp - 1._CUSTOM_REAL)) / &
                 (inv_tau_temp - alpha_only)
    eta_only = inv_tau_temp * (beta_only - alpha_only) / (inv_tau_temp - alpha_only)

    bar_A_1 = bar_A_0 * (gamma_only -alpha_only)
    bar_A_2 = 0._CUSTOM_REAL

    A_0 = bar_A_0
    A_1 = bar_A_1
    A_2 = bar_A_2

    A_zener = bar_A_0 * (tau_epsilon_temp * inv_tau_temp - 1._CUSTOM_REAL) * (eta_only - inv_tau_temp)

  else if (CPML_region_local == CPML_Z_ONLY_TEMP) then

    bar_A_0 = 1._CUSTOM_REAL / kappa_z
    alpha_only = beta_z
    beta_only = alpha_z

    gamma_only = (beta_only * (inv_tau_temp - alpha_only * tau_epsilon_temp * inv_tau_temp) + &
                  alpha_only * alpha_only * (tau_epsilon_temp * inv_tau_temp - 1._CUSTOM_REAL)) &
                  / (inv_tau_temp - alpha_only)
    eta_only = inv_tau_temp * (beta_only - alpha_only) / (inv_tau_temp - alpha_only)


    bar_A_1 = 0._CUSTOM_REAL
    bar_A_2 = bar_A_0* (gamma_only - alpha_only)

    A_0 = bar_A_0
    A_1 = bar_A_1
    A_2 = bar_A_2
    A_zener = bar_A_0 * (tau_epsilon_temp * inv_tau_temp - 1._CUSTOM_REAL) * (eta_only - inv_tau_temp)

  else

    call stop_the_code('error: CPML_region_local is undefined in lik_parameter_computation_viscoelastic()')

  endif

  end subroutine lik_parameter_computation_viscoelastic

!========================================================================

  subroutine l_parameter_computation(DT,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                     CPML_region_local,A_0,A_1,A_2,A_3,A_4, &
                                     bb_1,coef0_1,coef1_1,coef2_1,bb_2,coef0_2,coef1_2,coef2_2)

  use constants, only: CPML_X_ONLY,CPML_Z_ONLY,CPML_XZ,CUSTOM_REAL

  implicit none

  double precision, intent(in) :: DT
  double precision, intent(in) :: kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z
  integer, intent(in) :: CPML_region_local

  double precision, intent(out) :: A_0, A_1, A_2, A_3, A_4
  double precision, intent(out) :: bb_1,coef0_1,coef1_1,coef2_1,bb_2,coef0_2,coef1_2,coef2_2

  ! local variables
  double precision :: bar_A_0, bar_A_1, bar_A_2, bar_A_3, bar_A_4,gamma_x,gamma_z


  if (CPML_region_local == CPML_XZ) then
    bar_A_0 = kappa_x * kappa_z

    gamma_x = (alpha_x * alpha_z + alpha_x**2 + 2._CUSTOM_REAL * beta_x * beta_z - &
              2._CUSTOM_REAL * alpha_x * (beta_x + beta_z)) / (alpha_z - alpha_x)
    gamma_z = (alpha_x * alpha_z + alpha_z**2 + 2._CUSTOM_REAL * beta_x * beta_z - &
              2._CUSTOM_REAL * alpha_z * (beta_x + beta_z)) / (alpha_x - alpha_z)
    bar_A_1 = bar_A_0 /2._CUSTOM_REAL * (gamma_x - alpha_x + gamma_z - alpha_z)
    bar_A_2 = bar_A_0 /2._CUSTOM_REAL * (alpha_x**2 - gamma_x * alpha_x + &
              alpha_z**2 - gamma_z * alpha_z)

    A_0 = bar_A_0
    A_1 = bar_A_1
    A_2 = bar_A_2

    bar_A_3 = bar_A_0 /2._CUSTOM_REAL * alpha_x**2 * (gamma_x - alpha_x)
    bar_A_4 = bar_A_0 /2._CUSTOM_REAL * alpha_z**2 * (gamma_z - alpha_z)

    A_3 = bar_A_3
    A_4 = bar_A_4
  else if (CPML_region_local == CPML_X_ONLY) then
    bar_A_0 = kappa_x
    bar_A_1 = bar_A_0 * (beta_x - alpha_x)
    bar_A_2 = - bar_A_0 * alpha_x * (beta_x - alpha_x)

    A_0 = bar_A_0
    A_1 = bar_A_1
    A_2 = bar_A_2

    bar_A_3 = bar_A_0 * alpha_x**2 * (beta_x - alpha_x)
    bar_A_4 = 0.d0

    A_3 = bar_A_3
    A_4 = bar_A_4
  else if (CPML_region_local == CPML_Z_ONLY) then
    bar_A_0 = kappa_z
    bar_A_1 = bar_A_0 * (beta_z - alpha_z)
    bar_A_2 = - bar_A_0 * alpha_z * (beta_z - alpha_z)

    A_0 = bar_A_0
    A_1 = bar_A_1
    A_2 = bar_A_2

    bar_A_3 = 0.d0
    bar_A_4 = bar_A_0 * alpha_z**2 * (beta_z - alpha_z)

    A_3 = bar_A_3
    A_4 = bar_A_4
  endif

  bb_1 = alpha_x
  call compute_coef_convolution(bb_1,DT,coef0_1,coef1_1,coef2_1)

  bb_2 = alpha_z
  call compute_coef_convolution(bb_2,DT,coef0_2,coef1_2,coef2_2)

  end subroutine l_parameter_computation

!=====================================================================

  subroutine rebuild_value_on_PML_interface_acoustic(it,b_potential_acoustic,b_potential_dot_acoustic)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only: nspec,nglob,ispec_is_acoustic,any_acoustic,ibool,nglob_interface,point_interface

  ! PML arrays
  use specfem_par, only: ispec_is_PML,pml_interface_history_potential_dot,pml_interface_history_potential

  implicit none

  integer,intent(in) :: it

  real(kind=CUSTOM_REAL), dimension(nglob),intent(out) :: b_potential_acoustic,b_potential_dot_acoustic

  !local variables
  integer :: i,j,ispec

  if (any_acoustic) then
    do ispec = 1,nspec
      if (ispec_is_acoustic(ispec) .and. ispec_is_PML(ispec)) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            b_potential_dot_acoustic(ibool(i,j,ispec)) = 0._CUSTOM_REAL
            b_potential_acoustic(ibool(i,j,ispec)) = 0._CUSTOM_REAL
          enddo
        enddo
      endif
    enddo

    if (nglob_interface > 0) then
      do i = 1, nglob_interface
        b_potential_dot_acoustic(point_interface(i)) = pml_interface_history_potential_dot(i,it)
        b_potential_acoustic(point_interface(i)) = pml_interface_history_potential(i,it)
      enddo
    endif
  endif

  end subroutine rebuild_value_on_PML_interface_acoustic

!=====================================================================

  subroutine rebuild_value_on_PML_interface_acoustic_accel(it,b_potential_dot_dot_acoustic)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only: any_acoustic,nglob,nglob_interface,point_interface

  ! PML arrays
  use specfem_par, only: pml_interface_history_potential_dot_dot

  implicit none

  integer,intent(in) :: it

  real(kind=CUSTOM_REAL), dimension(nglob),intent(out) :: b_potential_dot_dot_acoustic

  !local variables
  integer :: i

  if (any_acoustic .and. nglob_interface > 0) then
    do i = 1, nglob_interface
      b_potential_dot_dot_acoustic(point_interface(i)) = pml_interface_history_potential_dot_dot(i,it)
    enddo
  endif

  end subroutine rebuild_value_on_PML_interface_acoustic_accel


!========================================================================

  subroutine rebuild_value_on_PML_interface_viscoelastic(it)

  use constants, only: NGLLX,NGLLZ,CUSTOM_REAL

  use specfem_par, only: nspec,ispec_is_elastic,any_elastic,ibool,nglob_interface,point_interface, &
                         b_veloc_elastic,b_displ_elastic
  ! PML arrays
  use specfem_par, only: ispec_is_PML,pml_interface_history_veloc,pml_interface_history_displ

  implicit none

  integer,intent(in) :: it

  !local variables
  integer :: i,j,ispec

  if (any_elastic) then
    do ispec = 1,nspec
      if (ispec_is_elastic(ispec) .and. ispec_is_PML(ispec)) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            b_veloc_elastic(:,ibool(i,j,ispec)) = 0._CUSTOM_REAL
            b_displ_elastic(:,ibool(i,j,ispec)) = 0._CUSTOM_REAL
          enddo
        enddo
      endif
    enddo

    if (nglob_interface > 0) then
      do i = 1, nglob_interface
        b_veloc_elastic(1,point_interface(i)) = pml_interface_history_veloc(1,i,it)
        b_veloc_elastic(2,point_interface(i)) = pml_interface_history_veloc(2,i,it)

        b_displ_elastic(1,point_interface(i)) = pml_interface_history_displ(1,i,it)
        b_displ_elastic(2,point_interface(i)) = pml_interface_history_displ(2,i,it)
      enddo
    endif
  endif

  end subroutine rebuild_value_on_PML_interface_viscoelastic

!========================================================================

  subroutine rebuild_value_on_PML_interface_viscoelastic_accel(it)

  use specfem_par, only: any_elastic,nglob_interface,point_interface, &
                         b_accel_elastic
  ! PML arrays
  use specfem_par, only: pml_interface_history_accel

  implicit none

  integer,intent(in) :: it

  !local variables
  integer :: i

  if (any_elastic .and. nglob_interface > 0) then
    do i = 1, nglob_interface
      b_accel_elastic(1,point_interface(i)) = pml_interface_history_accel(1,i,it)
      b_accel_elastic(2,point_interface(i)) = pml_interface_history_accel(2,i,it)
    enddo
  endif

  end subroutine rebuild_value_on_PML_interface_viscoelastic_accel

!========================================================================

  subroutine pml_boundary_acoustic(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                   potential_acoustic,potential_acoustic_old)

! The outer boundary condition to use for PML elements in fluid layers is Neumann for the potential
! because we need Dirichlet conditions for the displacement vector, which means Neumann for the potential.
! Thus, there is nothing to enforce explicitly here.
! There is something to enforce explicitly only in the case of elastic elements, for which a Dirichlet
! condition is needed for the displacement vector, which is the vectorial unknown for these elements.

!! DK DK this paragraph seems to be from Zhinan or from ChangHua:
! However, enforcing explicitly potential_dot_dot_acoustic, potential_dot_acoustic, potential_acoustic
! to be zero on outer boundary of PML help to improve the accuracy of absorbing low-frequency wave components
! in case of long-time simulation.

  use constants, only: CUSTOM_REAL

  use specfem_par, only: nglob_acoustic,anyabs,PML_abs_points_acoustic,PML_nglob_abs_acoustic

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic),intent(inout) :: potential_dot_dot_acoustic,potential_dot_acoustic, &
                                                                     potential_acoustic,potential_acoustic_old

  ! local parameters
  integer :: i,iglob

  ! checks if anything to do
  if (.not. anyabs) return
  if (PML_nglob_abs_acoustic == 0) return

  ! set Dirichelet boundary condition on outer boundary of CFS-PML
  do i = 1,PML_nglob_abs_acoustic
    iglob = PML_abs_points_acoustic(i)
    potential_acoustic_old(iglob) = 0._CUSTOM_REAL
    potential_acoustic(iglob) = 0._CUSTOM_REAL
    potential_dot_acoustic(iglob) = 0._CUSTOM_REAL
    potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
  enddo

  end subroutine pml_boundary_acoustic

!========================================================================

  subroutine pml_boundary_elastic(accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old)

  use constants, only: CUSTOM_REAL,NDIM

  use specfem_par, only: nglob_elastic,anyabs,PML_abs_points_elastic,PML_nglob_abs_elastic

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_elastic),intent(inout) :: accel_elastic,veloc_elastic,displ_elastic, &
                                                                         displ_elastic_old

  ! local parameters
  integer :: i,iglob

  ! checks if anything to do
  if (.not. anyabs) return
  if (PML_nglob_abs_elastic == 0) return

  !set Dirichlet boundary condition on outer boundary of PML

  ! we have to put Dirichlet on the boundary of the PML
  do i = 1,PML_nglob_abs_elastic
    iglob = PML_abs_points_elastic(i)
    displ_elastic_old(:,iglob) = 0._CUSTOM_REAL
    displ_elastic(:,iglob) = 0._CUSTOM_REAL
    veloc_elastic(:,iglob) = 0._CUSTOM_REAL
    accel_elastic(:,iglob) = 0._CUSTOM_REAL
  enddo

  end subroutine pml_boundary_elastic

!========================================================================

  subroutine decompose_rational_fraction_PML(N_PMLSF,beta,alpha)

  use constants, only: CUSTOM_REAL

  implicit none

  integer, intent(in) :: N_PMLSF
  double precision, dimension(N_PMLSF), intent(inout) :: beta,alpha

  !local parameter
  integer :: i
  double precision, dimension(N_PMLSF) ::gamma_de

  if (N_PMLSF < 0) then
    write(*,*)'the number of PML Stretching function should be greater than 0'
    call stop_the_code('error: stopping the code')
  endif


  if (N_PMLSF == 1) then
    ! there is nothing to do in this case
  else if (N_PMLSF == 2) then
    gamma_de(1) = (alpha(1) * alpha(2) + alpha(1)**2 + 2._CUSTOM_REAL * beta(1) * beta(2) - &
                  2._CUSTOM_REAL * alpha(1) * (beta(1) + beta(2))) / (alpha(2) - alpha(1))
    gamma_de(2) = (alpha(1) * alpha(2) + alpha(2)**2 + 2._CUSTOM_REAL * beta(1) * beta(2) - &
                  2._CUSTOM_REAL * alpha(2) * (beta(1) + beta(2))) / (alpha(1)-alpha(2))
  endif

  do i = 1,N_PMLSF
    beta(i) = gamma_de(i)
  enddo

 end subroutine decompose_rational_fraction_PML
 !=================================================
