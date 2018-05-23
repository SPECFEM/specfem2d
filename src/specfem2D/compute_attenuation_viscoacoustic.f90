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

! for viscoacoustic solver

!! DK DK QUENTIN visco begin

    ! ----------------------------------------------------------------------

  subroutine compute_attenuation_acoustic_integration(potential_acoustic,ispec_is_acoustic,PML_BOUNDARY_CONDITIONS,iphase,dot_e1)

  ! updates memory variable in viscoacoustic simulation

  ! compute forces for the elastic elements
  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,TWO,ALPHA_LDDRK,BETA_LDDRK,C_LDDRK,ZERO,ONE

  use specfem_par, only: nglob,nspec,N_SLS, &
                         ibool,xix,xiz,gammax,gammaz,hprime_xx,hprime_zz,ispec_is_PML, &
                         NGLJ, assign_external_model,AXISYM, &
                         nspec_inner_acoustic, nspec_outer_acoustic, &
                         kmato, phase_ispec_inner_acoustic,jacobian,rhoext,density, &
                         is_on_the_axis,hprimeBar_xx,hprimeBarwglj_xx,hprimewgll_zz,hprimewgll_xx, &
                         wxgll,wzgll,wxglj,xiglj,coord,iglob_is_forced,nglob_att

  implicit none

! update the memory variables using a convolution or using a differential equation
! (tests made by Ting Yu and also by Zhinan Xie, CNRS Marseille, France, show that it is better to leave it to .true.)
  logical, parameter :: CONVOLUTION_MEMORY_VARIABLES = .true.

  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: potential_acoustic

  logical,dimension(nspec),intent(in) :: ispec_is_acoustic

  ! CPML coefficients and memory variables
  logical,intent(in) :: PML_BOUNDARY_CONDITIONS

  integer :: iphase

  ! local variables
  integer :: ispec
  integer :: i,j

  real(kind=CUSTOM_REAL), dimension(nglob_att,N_SLS) :: dot_e1

  integer :: k,iglob

  ! spatial derivatives
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: dux_dxi,dux_dgamma
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: dux_dxl,dux_dzl

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: potential_elem
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2!,tempx3
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL), dimension(6,NGLLX,NGLLZ) :: deriv
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) :: rhol,fac
  real(kind=CUSTOM_REAL) :: temp1l,temp2l

  integer :: num_elements,ispec_p

  ! choses inner/outer elements
  if (iphase == 1) then
    num_elements = nspec_outer_acoustic
  else
    num_elements = nspec_inner_acoustic
  endif

  ! loop over relaxation mechanisms
  !do i_sls = 1,N_SLS

  ! If Newmark scheme, we only compute the displacement divergence that
  ! is the same for each SLS
  !if (time_stepping_scheme == 1 .and. i_sls > 1) cycle

  ! loop over spectral elements
  do ispec_p = 1,num_elements

    ! returns element id from stored element list
    ispec = phase_ispec_inner_acoustic(ispec_p,iphase)

    ! only for acoustic spectral elements
    if (.not. ispec_is_acoustic(ispec)) cycle
    if ((.not. PML_BOUNDARY_CONDITIONS) .or. (PML_BOUNDARY_CONDITIONS .and. (.not. ispec_is_PML(ispec)))) then

    ! gets local potential for element
    rhol = density(1,kmato(ispec))
    do j = 1,NGLLZ
      do i = 1,NGLLX

        ! convention to indicate that Q = 9999 i.e. that there is no viscoacousticity at that GLL point
        !if (inv_tau_sigma_nu1(i,j,ispec,i_sls) < 0.) cycle

        iglob = ibool(i,j,ispec)
        potential_elem(i,j) = potential_acoustic(iglob)

        ! stores local array for element xi/gamma/jacobian (for better performance)
        deriv(1,i,j) = xix(i,j,ispec)
        deriv(2,i,j) = xiz(i,j,ispec)
        deriv(3,i,j) = gammax(i,j,ispec)
        deriv(4,i,j) = gammaz(i,j,ispec)
        deriv(5,i,j) = jacobian(i,j,ispec)

        !! EQ:theta_n_u * phinu1 - e1(i,j,ispec,i_sls) * tauinvnu

        ! if external density model
        if (assign_external_model) then
          rhol = rhoext(i,j,ispec)
        endif
        deriv(6,i,j) = jacobian(i,j,ispec) / rhol

      enddo
    enddo

    ! first double loop over GLL points to compute and store gradients
    call mxm_2comp_singleA(dux_dxi,dux_dgamma,potential_elem,hprime_xx,hprime_zz)

    ! AXISYM case overwrites dux_dxi
    if (AXISYM) then
      if (is_on_the_axis(ispec)) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            ! derivative along x and along z
            dux_dxi(i,j) = 0._CUSTOM_REAL
            do k = 1,NGLLX
              dux_dxi(i,j) = dux_dxi(i,j) + potential_elem(k,j) * hprimeBar_xx(i,k)
            enddo
          enddo
        enddo
      endif
    endif

    ! gets derivatives of ux and uz with respect to x and z
    do j = 1,NGLLZ
      do i = 1,NGLLX
        xixl = deriv(1,i,j)
        xizl = deriv(2,i,j)
        gammaxl = deriv(3,i,j)
        gammazl = deriv(4,i,j)

        ! derivatives of potential
        dux_dxl(i,j) = dux_dxi(i,j) * xixl + dux_dgamma(i,j) * gammaxl
        dux_dzl(i,j) = dux_dxi(i,j) * xizl + dux_dgamma(i,j) * gammazl
      enddo
    enddo

    ! AXISYM case overwrite dux_dxl
    if (AXISYM) then
      if (is_on_the_axis(ispec)) then
        ! dchi/dr=rho * u_r=0 on the axis
        ! i == 1
        do j = 1,NGLLZ
          dux_dxl(1,j) = 0._CUSTOM_REAL
        enddo
      endif
    endif

    ! first double loop to compute gradient
    if (AXISYM) then
      ! AXISYM case
      if (is_on_the_axis(ispec)) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            xixl = deriv(1,i,j)
            xizl = deriv(2,i,j)
            gammaxl = deriv(3,i,j)
            gammazl = deriv(4,i,j)
            jacobianl = deriv(5,i,j)
            fac = deriv(6,i,j) ! jacobian/rho

            if (i == 1) then
              ! dchi/dr=rho * u_r=0 on the axis
              dux_dxl(i,j) = 0._CUSTOM_REAL
              r_xiplus1(i,j) = gammaz(i,j,ispec) * jacobianl
            else
              r_xiplus1(i,j) = coord(1,ibool(i,j,ispec))/(xiglj(i) + ONE)
            endif
            tempx1(i,j) = r_xiplus1(i,j) * fac * (xixl * dux_dxl(i,j) + xizl * dux_dzl(i,j))
            tempx2(i,j) = r_xiplus1(i,j) * fac * (gammaxl * dux_dxl(i,j) + gammazl * dux_dzl(i,j))
          enddo
        enddo
      else
        do j = 1,NGLLZ
          do i = 1,NGLLX
            xixl = deriv(1,i,j)
            xizl = deriv(2,i,j)
            gammaxl = deriv(3,i,j)
            gammazl = deriv(4,i,j)
            jacobianl = deriv(5,i,j)
            fac = deriv(6,i,j) ! jacobian/rho

            tempx1(i,j) = coord(1,ibool(i,j,ispec)) * fac * (xixl * dux_dxl(i,j) + xizl * dux_dzl(i,j))
            tempx2(i,j) = coord(1,ibool(i,j,ispec)) * fac * (gammaxl * dux_dxl(i,j) + gammazl * dux_dzl(i,j))
          enddo
        enddo
      endif
    else
      ! default case
      do j = 1,NGLLZ
        do i = 1,NGLLX
          xixl = deriv(1,i,j)
          xizl = deriv(2,i,j)
          gammaxl = deriv(3,i,j)
          gammazl = deriv(4,i,j)
          jacobianl = deriv(5,i,j)
          fac = deriv(6,i,j) ! jacobian/rho

          tempx1(i,j) = fac * (xixl * dux_dxl(i,j) + xizl * dux_dzl(i,j))
          tempx2(i,j) = fac * (gammaxl * dux_dxl(i,j) + gammazl * dux_dzl(i,j))
        enddo
      enddo
    endif

!
! second double-loop over GLL to compute all the terms
!

    ! along x direction and z direction
    ! and assemble the contributions
    if (AXISYM) then
      ! axisymmetric case
      if (is_on_the_axis(ispec)) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            if (.not. iglob_is_forced(iglob)) then
              ! assembles the contributions
              temp1l = 0._CUSTOM_REAL
              temp2l = 0._CUSTOM_REAL
              do k = 1,NGLLX
                temp1l = temp1l + tempx1(k,j) * hprimeBarwglj_xx(k,i)
                temp2l = temp2l + tempx2(i,k) * hprimewgll_zz(k,j)
              enddo
              ! sums contributions from each element to the global values
              dot_e1(iglob,:) = dot_e1(iglob,:) &
                                                  - (wzgll(j) * temp1l + wxglj(i) * temp2l)
            endif
          enddo
        enddo
      else
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            if (.not. iglob_is_forced(iglob)) then
              ! assembles the contributions
              temp1l = 0._CUSTOM_REAL
              temp2l = 0._CUSTOM_REAL
              do k = 1,NGLLX
                temp1l = temp1l + tempx1(k,j) * hprimewgll_xx(k,i)
                temp2l = temp2l + tempx2(i,k) * hprimewgll_zz(k,j)
              enddo
              ! sums contributions from each element to the global values
              dot_e1(iglob,:) = dot_e1(iglob,:) &
                                                  - (wzgll(j) * temp1l + wxgll(i) * temp2l)
            endif
          enddo
        enddo
      endif
    else
      ! default case
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          if (.not. iglob_is_forced(iglob)) then
            ! assembles the contributions
            temp1l = 0._CUSTOM_REAL
            temp2l = 0._CUSTOM_REAL
            do k = 1,NGLLX
              temp1l = temp1l + tempx1(k,j) * hprimewgll_xx(k,i)
              temp2l = temp2l + tempx2(i,k) * hprimewgll_zz(k,j)
            enddo
            ! sums contributions from each element to the global values
            dot_e1(iglob,:) = dot_e1(iglob,:) &
                                                - (wzgll(j) * temp1l + wxgll(i) * temp2l)
          endif
        enddo
      enddo
    endif

   endif

  enddo ! end of loop over all spectral elements

  !enddo ! end of loop over all relaxation mechanisms

  contains

!---------------------------------------------------------------------------------------

  subroutine mxm_2comp_singleA(x,z,A,B,C)

! matrix x matrix multiplication, merging 2 loops for x = A^t B^t and z = A C^t
!
! index notation:
! general matrix multiplication: uij = (A B)ij = Aik Bkj
!                          here: xij = (A^t B^t)ij = Akj Bik = (B A)ij
!                                zij = (A C^t)ij = Aik Cjk
!
! original loops:
!
!      do j = 1,NGLLZ
!        do i = 1,NGLLX
!          ! derivative along x and along z
!          dux_dxi(i,j) = 0._CUSTOM_REAL
!          dux_dgamma(i,j) = 0._CUSTOM_REAL
!
!          ! first double loop over GLL points to compute and store gradients
!          ! we can merge the two loops because NGLLX == NGLLZ
!          do k = 1,NGLLX
!            dux_dxi(i,j) = dux_dxi(i,j) + potential_elem(k,j) * hprime_xx(i,k)
!            dux_dgamma(i,j) = dux_dgamma(i,j) + potential_elem(i,k) * hprime_zz(j,k)
!          enddo
!        enddo
!      enddo

  use constants, only: NGLLX,NGLLZ,CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(out) :: x,z
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(in) :: A,B,C

  ! local parameters
  integer :: i,j,k

  select case(NGLLX)
  case (5)
    do j = 1,5
      do i = 1,5
        ! loop unrolling
        x(i,j) = A(1,j) * B(i,1) + A(2,j) * B(i,2) + A(3,j) * B(i,3) + A(4,j) * B(i,4) + A(5,j) * B(i,5)
        z(i,j) = A(i,1) * C(j,1) + A(i,2) * C(j,2) + A(i,3) * C(j,3) + A(i,4) * C(j,4) + A(i,5) * C(j,5)
      enddo
    enddo

  case default
    do j = 1,NGLLZ
      do i = 1,NGLLX
        x(i,j) = 0._CUSTOM_REAL
        z(i,j) = 0._CUSTOM_REAL
        do k = 1,NGLLX
          x(i,j) = x(i,j) + A(k,j) * B(i,k)
          z(i,j) = z(i,j) + A(i,k) * C(j,k)
        enddo
      enddo
    enddo
  end select

  end subroutine mxm_2comp_singleA

  end subroutine compute_attenuation_acoustic_integration

  ! -----------------------------------------------------------------------------------------------------------------------

  subroutine update_memory_var_acous_weak_form(dot_e1)

  ! compute forces for the elastic elements
  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,TWO,ALPHA_LDDRK,BETA_LDDRK,C_LDDRK,ZERO

  use specfem_par, only: N_SLS, &
                         time_stepping_scheme,i_stage,deltat, &
                         e1_LDDRK_acous, e1_initial_rk_acous, e1_force_RK_acous, &
                         e1_acous,e1_acous_temp,initialfield,it, &
                         A_newmark_e1, B_newmark_e1, phi_nu1, inv_tau_sigma_nu1, nspec, &
                         PML_BOUNDARY_CONDITIONS, ispec_is_PML, ispec_is_acoustic, ibool, &
                         dot_e1_old,nglob_acoustic

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic,N_SLS),intent(inout) :: dot_e1

  ! Local variables
  integer :: i_sls, iglob, i, j, ispec
  real(kind=CUSTOM_REAL) :: weight_rk, phinu1, tauinvnu1, temp, coef1

  ! loop over relaxation mechanisms
  do i_sls = 1,N_SLS

  ! Time update
  select case (time_stepping_scheme)

        case (1)

          ! Newmark
          ! Set newmark coefficients
          if (it == 1 .and. i_stage == 1) then

          ! loop over spectral elements
          do ispec = 1,nspec

             if (.not. ispec_is_acoustic(ispec)) cycle

             if ((.not. PML_BOUNDARY_CONDITIONS) .or. (PML_BOUNDARY_CONDITIONS &
                        .and. (.not. ispec_is_PML(ispec)))) then

                  do j = 1,NGLLZ
                    do i = 1,NGLLX

                    iglob = ibool(i,j,ispec)

                    phinu1    = phi_nu1(i,j,ispec,i_sls)
                    tauinvnu1 = inv_tau_sigma_nu1(i,j,ispec,i_sls)
                    temp      = exp(- 0.5d0 * tauinvnu1 * deltat)
                    coef1     = (1.d0 - temp) / tauinvnu1

                    A_newmark_e1(iglob,i_sls) = temp
                    B_newmark_e1(iglob,i_sls) = phinu1 * coef1

                    enddo
                  enddo

             endif

          enddo

          endif

! update the memory variables using a convolution or using a differential equation
! From Zhinan Xie and Dimitri Komatitsch:
! For cases in which a value of tau_sigma is small, then its inverse is large,
! which may result in a in stiff ordinary differential equation to solve;
! in such a case, resorting to the convolution formulation is better.
!! DK DK inlined this for speed            call compute_coef_convolution(tauinvnu1,deltat,coef0,coef1,coef2)
            !temp = exp(- 0.5d0 * tauinvnu1 * deltat)
            !coef1 = (1.d0 - temp) / tauinvnu1

            e1_acous(:,i_sls) = (A_newmark_e1(:,i_sls)**2) * e1_acous(:,i_sls) &
                + B_newmark_e1(:,i_sls) * (dot_e1(:,1) + A_newmark_e1(:,i_sls) * dot_e1_old(:,1))

            dot_e1_old = dot_e1

        case (2)
          ! LDDRK
          ! update e1, e11, e13 in ADE formation with fourth-order LDDRK scheme
          e1_LDDRK_acous(:,i_sls) = ALPHA_LDDRK(i_stage) * e1_LDDRK_acous(:,i_sls) + deltat * dot_e1(:,i_sls)

          if (i_stage == 1 .and. it == 1 .and. (.not. initialfield) .and. .false.) then
            !! DK DK this should be vectorized
            e1_acous_temp(:,i_sls) = e1_acous_temp(:,i_sls) + &
                                     BETA_LDDRK(i_stage) * e1_LDDRK_acous(:,i_sls)
            e1_acous(:,i_sls) = e1_acous_temp(:,i_sls)
          else
            e1_acous(:,i_sls) = e1_acous(:,i_sls) + BETA_LDDRK(i_stage) * e1_LDDRK_acous(:,i_sls)
          endif

        case (3)
          ! Runge-Kutta
          ! update e1, e11, e13 in ADE formation with classical fourth-order Runge-Kutta scheme
          e1_force_RK_acous(:,i_sls,i_stage) = deltat * dot_e1(:,i_sls)

          if (i_stage == 1 .or. i_stage == 2 .or. i_stage == 3) then
            if (i_stage == 1) weight_rk = 0.5_CUSTOM_REAL
            if (i_stage == 2) weight_rk = 0.5_CUSTOM_REAL
            if (i_stage == 3) weight_rk = 1._CUSTOM_REAL

            if (i_stage == 1) e1_initial_rk_acous(:,i_sls) = e1_acous(:,i_sls)
            e1_acous(:,i_sls) = e1_initial_rk_acous(:,i_sls) &
                + weight_rk * e1_force_RK_acous(:,i_sls,i_stage)
          else if (i_stage == 4) then
            e1_acous(:,i_sls) = e1_initial_rk_acous(:,i_sls) + 1._CUSTOM_REAL / 6._CUSTOM_REAL * &
                                  (e1_force_RK_acous(:,i_sls,1) + 2._CUSTOM_REAL * e1_force_RK_acous(:,i_sls,2) + &
                                   2._CUSTOM_REAL * e1_force_RK_acous(:,i_sls,3) + e1_force_RK_acous(:,i_sls,4))
          endif

        case default
          call stop_the_code('Time stepping scheme not implemented yet in viscoacoustic attenuation update')

  end select

  enddo ! end loop over all relaxation mechanisms

  dot_e1(:,:) = ZERO

  end subroutine update_memory_var_acous_weak_form

  ! -----------------------------------------------------------------------------------------------------------------------

  subroutine get_attenuation_forces_strong_form(sum_forces,sum_forces_old,forces_attenuation,i,j,ispec,iglob,e1_acous_sf)

  ! updates memory variable in viscoacoustic simulation
  ! and get the attenuation contribution

  ! compute forces for the elastic elements
  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,TWO,ALPHA_LDDRK,BETA_LDDRK,C_LDDRK

  use specfem_par, only: PML_BOUNDARY_CONDITIONS,N_SLS,nspec_ATT_ac, &
                         ispec_is_PML, &
                         phi_nu1, inv_tau_sigma_nu1,time_stepping_scheme,i_stage,deltat, &
                         e1_acous, e1_LDDRK_acous, e1_initial_rk_acous, e1_force_RK_acous, &
                         A_newmark_e1_sf, B_newmark_e1_sf

  implicit none


  real(kind=CUSTOM_REAL), intent(in)    :: sum_forces
  real(kind=CUSTOM_REAL), intent(inout) :: sum_forces_old
  real(kind=CUSTOM_REAL), intent(out)   :: forces_attenuation
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLZ,nspec_ATT_ac), intent(inout):: e1_acous_sf
  integer :: i,j,ispec,iglob

  ! local variables

! update the memory variables using a convolution or using a differential equation
! (tests made by Ting Yu and also by Zhinan Xie, CNRS Marseille, France, show
! that it is better to leave it to .true.)
  logical, parameter :: CONVOLUTION_MEMORY_VARIABLES = .true.

  ! for attenuation
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: phinu1,a_newmark
  double precision, dimension(N_SLS) :: tauinvnu1

  ! temporary RK4 variable
  real(kind=CUSTOM_REAL) :: weight_rk

  if (.not. CONVOLUTION_MEMORY_VARIABLES) &
    call stop_the_code('CONVOLUTION_MEMORY_VARIABLES == .false. is not accurate enough and has been discontinued for now')

  forces_attenuation = 0._CUSTOM_REAL

  if ( PML_BOUNDARY_CONDITIONS .and. ispec_is_PML(ispec)) return

  ! convention to indicate that Q = 9999 i.e. that there is no viscoacousticity at that GLL point
  if (inv_tau_sigma_nu1(i,j,ispec,1) < 0.) return

  if (time_stepping_scheme /= 1) then
    phinu1(:)    = phi_nu1(i,j,ispec,:)
    tauinvnu1(:) = inv_tau_sigma_nu1(i,j,ispec,:)
  endif

  ! update e1 in convolution formulation with modified recursive convolution scheme on basis of
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
    a_newmark = A_newmark_e1_sf(:,i,j,ispec)
    e1_acous_sf(:,i,j,ispec) = a_newmark * a_newmark * e1_acous_sf(:,i,j,ispec) + &
                               B_newmark_e1_sf(:,i,j,ispec) * (sum_forces + a_newmark * sum_forces_old)
    forces_attenuation = sum(e1_acous_sf(:,i,j,ispec))
    sum_forces_old = sum_forces
!         else
!           stop 'CONVOLUTION_MEMORY_VARIABLES == .false. is not accurate enough and has been discontinued for now'
!           e1(i,j,ispec,i_sls) = e1(i,j,ispec,i_sls) + deltat * (- e1(i,j,ispec,i_sls)*tauinvnu1 + phinu1 * theta_n_u)
!         endif

  case (2)
    ! LDDRK
    ! update e1, e11, e13 in ADE formation with fourth-order LDDRK scheme
    e1_LDDRK_acous(iglob,:) = ALPHA_LDDRK(i_stage) * e1_LDDRK_acous(iglob,:) + &
                                  deltat * (sum_forces * phinu1 - e1_acous(iglob,:) * tauinvnu1)
    e1_acous(iglob,:) = e1_acous(iglob,:) + BETA_LDDRK(i_stage) * e1_LDDRK_acous(iglob,:)

    forces_attenuation = sum(e1_acous(iglob,:))

  case (3)
    ! Runge-Kutta
    ! update e1, e11, e13 in ADE formation with classical fourth-order Runge-Kutta scheme
    e1_force_RK_acous(iglob,:,i_stage) = deltat * (sum_forces * phinu1 - e1_acous(iglob,:) * tauinvnu1)

    if (i_stage == 1 .or. i_stage == 2 .or. i_stage == 3) then
      if (i_stage == 1) weight_rk = 0.5_CUSTOM_REAL
      if (i_stage == 2) weight_rk = 0.5_CUSTOM_REAL
      if (i_stage == 3) weight_rk = 1._CUSTOM_REAL

      if (i_stage == 1) e1_initial_rk_acous(iglob,:) = e1_acous(iglob,:)
      e1_acous(iglob,:) = e1_initial_rk_acous(iglob,:) &
                              + weight_rk * e1_force_RK_acous(iglob,:,i_stage)
    else if (i_stage == 4) then
      e1_acous(iglob,:) = e1_initial_rk_acous(iglob,:) + 1._CUSTOM_REAL / 6._CUSTOM_REAL * &
                             (e1_force_RK_acous(iglob,:,1) + 2._CUSTOM_REAL * e1_force_RK_acous(iglob,:,2) + &
                              2._CUSTOM_REAL * e1_force_RK_acous(iglob,:,3) + e1_force_RK_acous(iglob,:,4))
    endif

    forces_attenuation = sum(e1_acous(iglob,:))

  case default
    call stop_the_code('Time stepping scheme not implemented yet in viscoacoustic attenuation update')
  end select

  end subroutine get_attenuation_forces_strong_form
