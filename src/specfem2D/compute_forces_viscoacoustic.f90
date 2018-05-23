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

  subroutine compute_forces_viscoacoustic(potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic, &
                                          PML_BOUNDARY_CONDITIONS,potential_acoustic_old,iphase,e1_acous_sf,sum_forces_old)

! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,CPML_X_ONLY,CPML_Z_ONLY,IRIGHT,ILEFT,IBOTTOM,ITOP, &
    ZERO,ONE,TWO,TWO_THIRDS, &
    ALPHA_LDDRK,BETA_LDDRK,C_LDDRK,USE_A_STRONG_FORMULATION_FOR_E1

  use specfem_par, only: nglob,nspec_ATT_ac, &
                         assign_external_model,ibool,kmato,ispec_is_acoustic, &
                         density,rhoext, &
                         xix,xiz,gammax,gammaz,jacobian, &
                         hprime_xx,hprimewgll_xx, &
                         hprime_zz,hprimewgll_zz,wxgll,wzgll, &
                         AXISYM,is_on_the_axis,coord,hprimeBar_xx,hprimeBarwglj_xx,xiglj,wxglj,ATTENUATION_VISCOACOUSTIC, &
                         N_SLS, iglob_is_forced,time_stepping_scheme,phi_nu1,inv_tau_sigma_nu1, &
                         e1_acous,dot_e1

  ! overlapping communication
  use specfem_par, only: nspec_inner_acoustic,nspec_outer_acoustic,phase_ispec_inner_acoustic

  ! PML arrays
  use specfem_par, only: ispec_is_PML

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob),intent(inout) :: potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: potential_dot_acoustic,potential_acoustic

  logical,intent(in) :: PML_BOUNDARY_CONDITIONS
  real(kind=CUSTOM_REAL), dimension(nglob) :: potential_acoustic_old

  real(kind=CUSTOM_REAL),dimension(N_SLS,NGLLX,NGLLZ,nspec_ATT_ac) :: e1_acous_sf
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ,nspec_ATT_ac) :: sum_forces_old

  integer,intent(in) :: iphase

  ! local parameters
  integer :: ispec,i,j,k,iglob

  ! spatial derivatives
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: dux_dxi,dux_dgamma
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: dux_dxl,dux_dzl

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: potential_elem
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2,tempx3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,N_SLS) :: tempx3_e1
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL), dimension(6,NGLLX,NGLLZ) :: deriv
  real(kind=CUSTOM_REAL), dimension(2,NGLLX,NGLLZ,N_SLS) :: deriv_e1
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) :: rhol,fac
  real(kind=CUSTOM_REAL) :: temp1l,temp2l,sum_forces,forces_attenuation

  ! local PML parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: potential_dot_dot_acoustic_PML

  integer :: num_elements,ispec_p

  integer :: i_sls

  ! choses inner/outer elements
  if (iphase == 1) then
    num_elements = nspec_outer_acoustic
  else
    num_elements = nspec_inner_acoustic
  endif

  ! loop over spectral elements
  do ispec_p = 1,num_elements

    ! returns element id from stored element list
    ispec = phase_ispec_inner_acoustic(ispec_p,iphase)

    ! only for acoustic spectral elements
    if (.not. ispec_is_acoustic(ispec)) cycle

    ! gets local potential for element
    rhol = density(1,kmato(ispec))
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        potential_elem(i,j) = potential_acoustic(iglob)

        ! stores local array for element xi/gamma/jacobian (for better performance)
        deriv(1,i,j) = xix(i,j,ispec)
        deriv(2,i,j) = xiz(i,j,ispec)
        deriv(3,i,j) = gammax(i,j,ispec)
        deriv(4,i,j) = gammaz(i,j,ispec)
        deriv(5,i,j) = jacobian(i,j,ispec)
        ! if external density model
        if (assign_external_model) then
          rhol = rhoext(i,j,ispec)
        endif
        deriv(6,i,j) = jacobian(i,j,ispec) / rhol

        if (ATTENUATION_VISCOACOUSTIC .and. (.not. USE_A_STRONG_FORMULATION_FOR_E1) .and. time_stepping_scheme > 1) then
                deriv_e1(1,i,j,:) = phi_nu1(i,j,ispec,:)
                deriv_e1(2,i,j,:) = inv_tau_sigma_nu1(i,j,ispec,:)
        endif

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

    ! derivative along x and along zbb
    if (PML_BOUNDARY_CONDITIONS) then
      call pml_compute_memory_variables_acoustic(ispec,nglob,potential_acoustic_old,dux_dxl,dux_dzl)
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

    ! first double loop over GLL points to compute and store gradients
    if (PML_BOUNDARY_CONDITIONS) then
      ! calculates contribution from each C-PML element to update acceleration
      call pml_compute_accel_contribution_acoustic(ispec,nglob, &
                                                   potential_acoustic,potential_acoustic_old,potential_dot_acoustic, &
                                                   potential_dot_dot_acoustic_PML,r_xiplus1)
    endif

!! DK DK QUENTIN visco begin
    if (ATTENUATION_VISCOACOUSTIC .and. (.not. USE_A_STRONG_FORMULATION_FOR_E1)) then
      ! attenuation is implemented following the memory variable formulation of
      ! Carcione et al., Wave propagation simulation in a linear viscoacoustic medium,
      ! Geophysical Journal, vol. 93, p. 393-407 (1988)

      tempx3    = 0._CUSTOM_REAL
      tempx3_e1 = 0._CUSTOM_REAL
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          ! loop on all the standard linear solids
          do i_sls = 1,N_SLS
            tempx3(i,j) = tempx3(i,j) + e1_acous(iglob,i_sls)
          enddo
          if (time_stepping_scheme > 1) tempx3_e1(i,j,:) = deriv(5,i,j)*(deriv_e1(2,i,j,:)/deriv_e1(1,i,j,:))*e1_acous(iglob,:)
          tempx3(i,j) = jacobian(i,j,ispec)  * tempx3(i,j)
        enddo
      enddo

    endif
!! DK DK QUENTIN visco end

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
              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                                  - (wzgll(j) * temp1l + wxglj(i) * temp2l)
              if (ATTENUATION_VISCOACOUSTIC .and. (.not. USE_A_STRONG_FORMULATION_FOR_E1) &
                   .and. time_stepping_scheme > 1) dot_e1(iglob,:) = dot_e1(iglob,:) &
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
              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                                  - (wzgll(j) * temp1l + wxgll(i) * temp2l)
              if (ATTENUATION_VISCOACOUSTIC .and. (.not. USE_A_STRONG_FORMULATION_FOR_E1) &
                  .and. time_stepping_scheme > 1) dot_e1(iglob,:) = dot_e1(iglob,:) &
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
            sum_forces = wzgll(j) * temp1l + wxgll(i) * temp2l
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - sum_forces
            if (ATTENUATION_VISCOACOUSTIC .and. (.not. USE_A_STRONG_FORMULATION_FOR_E1) &
                .and. time_stepping_scheme > 1) dot_e1(iglob,:) = dot_e1(iglob,:) - sum_forces
            if (ATTENUATION_VISCOACOUSTIC .and. USE_A_STRONG_FORMULATION_FOR_E1) then
              call get_attenuation_forces_strong_form(sum_forces,sum_forces_old(i,j,ispec), &
                                                      forces_attenuation,i,j,ispec,iglob,e1_acous_sf)
              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - forces_attenuation
            endif
          endif
        enddo
      enddo
    endif

    ! PML contribution
    if (PML_BOUNDARY_CONDITIONS) then
      if (ispec_is_PML(ispec)) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            if (.not. iglob_is_forced(iglob)) then
              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_dot_acoustic_PML(i,j)
            endif
          enddo
        enddo
      endif
    endif

!! DK DK QUENTIN visco begin
    if (ATTENUATION_VISCOACOUSTIC .and. (.not. USE_A_STRONG_FORMULATION_FOR_E1)) then
      if (AXISYM) then
        ! axisymmetric case
        if (is_on_the_axis(ispec)) then
          do j = 1,NGLLZ
            do i = 1,NGLLX
              ! sums contributions from each element to the global values
              iglob = ibool(i,j,ispec)
              if (.not. iglob_is_forced(iglob)) then
                potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + wzgll(j) * wxglj(i) * tempx3(i,j)
                ! loop over all relaxation mechanisms
                if (time_stepping_scheme > 1) dot_e1(iglob,:) = dot_e1(iglob,:) - wzgll(j) * wxglj(i) * tempx3_e1(i,j,:)
              endif
            enddo
          enddo
        else
          do j = 1,NGLLZ
            do i = 1,NGLLX
              ! sums contributions from each element to the global values
              iglob = ibool(i,j,ispec)
              if (.not. iglob_is_forced(iglob)) then
                potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + wzgll(j) * wxgll(i) * tempx3(i,j)
                ! loop over all relaxation mechanisms
                if (time_stepping_scheme > 1) dot_e1(iglob,:) = dot_e1(iglob,:) - wzgll(j) * wxgll(i) * tempx3_e1(i,j,:)
              endif
            enddo
          enddo
        endif
      else
        do j = 1,NGLLZ
          do i = 1,NGLLX
            ! sums contributions from each element to the global values
            iglob = ibool(i,j,ispec)
            if (.not. iglob_is_forced(iglob)) then
              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + wzgll(j) * wxgll(i) * tempx3(i,j)
              ! loop over all relaxation mechanisms
              if (time_stepping_scheme > 1) dot_e1(iglob,:) = dot_e1(iglob,:) - wzgll(j) * wxgll(i) * tempx3_e1(i,j,:)
            endif
          enddo
        enddo
      endif
    endif
!! DK DK QUENTIN visco end

  enddo ! end of loop over all spectral elements

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

  end subroutine compute_forces_viscoacoustic
