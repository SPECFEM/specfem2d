!========================================================================
!
!                            S P E C F E M 2 D
!                            -----------------
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

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJ, DO_LOOP_IJ, ENDDO_LOOP_IJ defined in config.fh
#include "config.fh"


  subroutine compute_forces_viscoacoustic(potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic, &
                                          PML_BOUNDARY_CONDITIONS,potential_acoustic_old,iphase,e1_acous_sf,sum_forces_old)

! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,CPML_X_ONLY,CPML_Z_ONLY,IRIGHT,ILEFT,IBOTTOM,ITOP, &
    ZERO,ONE,TWO,TWO_THIRDS, &
    USE_A_STRONG_FORMULATION_FOR_E1

  use specfem_par, only: nglob,nspec_ATT_ac, &
                         ibool,ispec_is_acoustic, &
                         deriv_mapping, &
                         hprime_xx,hprimewgll_xx, &
                         hprime_zz,hprimewgll_zz,wxgll,wzgll, &
                         AXISYM,is_on_the_axis,coord,hprimeBar_xx,hprimeBarwglj_xx,xiglj,wxglj,ATTENUATION_VISCOACOUSTIC, &
                         N_SLS, iglob_is_forced,time_stepping_scheme,phi_nu1,inv_tau_sigma_nu1, &
                         e1_acous,dot_e1

  ! overlapping communication
  use specfem_par, only: nspec_inner_acoustic,nspec_outer_acoustic,phase_ispec_inner_acoustic

  ! PML arrays
  use specfem_par, only: ispec_is_PML

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLSQUARE
#endif

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
#ifdef FORCE_VECTORIZATION
  integer :: ij
#endif

  ! spatial derivatives
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: dux_dxi,dux_dgamma
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: dux_dxl,dux_dzl

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: potential_elem
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2,tempx3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,N_SLS) :: tempx3_e1
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL), dimension(2,NGLLX,NGLLZ,N_SLS) :: deriv_e1
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) :: fac
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
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        potential_elem(i,j) = potential_acoustic(iglob)

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
            do k = 1,NGLJ
              dux_dxi(i,j) = dux_dxi(i,j) + potential_elem(k,j) * hprimeBar_xx(i,k)
            enddo
          enddo
        enddo
      endif
    endif

    ! gets derivatives of ux and uz with respect to x and z
    DO_LOOP_IJ
        xixl = deriv_mapping(1,INDEX_IJ,ispec)
        xizl = deriv_mapping(2,INDEX_IJ,ispec)
        gammaxl = deriv_mapping(3,INDEX_IJ,ispec)
        gammazl = deriv_mapping(4,INDEX_IJ,ispec)

        ! derivatives of potential
        dux_dxl(INDEX_IJ) = dux_dxi(INDEX_IJ) * xixl + dux_dgamma(INDEX_IJ) * gammaxl
        dux_dzl(INDEX_IJ) = dux_dxi(INDEX_IJ) * xizl + dux_dgamma(INDEX_IJ) * gammazl
    ENDDO_LOOP_IJ

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
            xixl = deriv_mapping(1,i,j,ispec)
            xizl = deriv_mapping(2,i,j,ispec)
            gammaxl = deriv_mapping(3,i,j,ispec)
            gammazl = deriv_mapping(4,i,j,ispec)
            jacobianl = deriv_mapping(5,i,j,ispec)
            fac = deriv_mapping(6,i,j,ispec) ! jacobian/rho

            if (i == 1) then
              ! dchi/dr=rho * u_r=0 on the axis
              dux_dxl(i,j) = 0._CUSTOM_REAL
              r_xiplus1(i,j) = gammazl * jacobianl
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
            xixl = deriv_mapping(1,i,j,ispec)
            xizl = deriv_mapping(2,i,j,ispec)
            gammaxl = deriv_mapping(3,i,j,ispec)
            gammazl = deriv_mapping(4,i,j,ispec)
            jacobianl = deriv_mapping(5,i,j,ispec)
            fac = deriv_mapping(6,i,j,ispec) ! jacobian/rho

            iglob = ibool(i,j,ispec)
            tempx1(i,j) = coord(1,iglob) * fac * (xixl * dux_dxl(i,j) + xizl * dux_dzl(i,j))
            tempx2(i,j) = coord(1,iglob) * fac * (gammaxl * dux_dxl(i,j) + gammazl * dux_dzl(i,j))
          enddo
        enddo
      endif
    else
      ! default case
      DO_LOOP_IJ
          xixl = deriv_mapping(1,INDEX_IJ,ispec)
          xizl = deriv_mapping(2,INDEX_IJ,ispec)
          gammaxl = deriv_mapping(3,INDEX_IJ,ispec)
          gammazl = deriv_mapping(4,INDEX_IJ,ispec)
          jacobianl = deriv_mapping(5,INDEX_IJ,ispec)
          fac = deriv_mapping(6,INDEX_IJ,ispec) ! jacobian/rho

          tempx1(INDEX_IJ) = fac * (xixl * dux_dxl(INDEX_IJ) + xizl * dux_dzl(INDEX_IJ))
          tempx2(INDEX_IJ) = fac * (gammaxl * dux_dxl(INDEX_IJ) + gammazl * dux_dzl(INDEX_IJ))
      ENDDO_LOOP_IJ
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
          jacobianl = deriv_mapping(5,i,j,ispec)

          ! loop on all the standard linear solids
          do i_sls = 1,N_SLS
            tempx3(i,j) = tempx3(i,j) + e1_acous(iglob,i_sls)
          enddo

          if (time_stepping_scheme > 1) then
            tempx3_e1(i,j,:) = jacobianl * (deriv_e1(2,i,j,:)/deriv_e1(1,i,j,:)) * e1_acous(iglob,:)
          endif

          tempx3(i,j) = jacobianl * tempx3(i,j)
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
              do k = 1,NGLJ
                temp1l = temp1l + tempx1(k,j) * hprimeBarwglj_xx(k,i)
                temp2l = temp2l + tempx2(i,k) * hprimewgll_zz(k,j)
              enddo
              ! sums contributions from each element to the global values
              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                                  - (wzgll(j) * temp1l + wxglj(i) * temp2l)
              if (ATTENUATION_VISCOACOUSTIC .and. (.not. USE_A_STRONG_FORMULATION_FOR_E1) .and. time_stepping_scheme > 1) &
                dot_e1(iglob,:) = dot_e1(iglob,:) - (wzgll(j) * temp1l + wxglj(i) * temp2l)
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
              if (ATTENUATION_VISCOACOUSTIC .and. (.not. USE_A_STRONG_FORMULATION_FOR_E1) .and. time_stepping_scheme > 1) &
                dot_e1(iglob,:) = dot_e1(iglob,:) - (wzgll(j) * temp1l + wxgll(i) * temp2l)
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

            if (ATTENUATION_VISCOACOUSTIC .and. (.not. USE_A_STRONG_FORMULATION_FOR_E1) .and. time_stepping_scheme > 1) &
              dot_e1(iglob,:) = dot_e1(iglob,:) - sum_forces

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
