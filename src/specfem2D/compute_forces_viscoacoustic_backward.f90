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

  subroutine compute_forces_viscoacoustic_backward(b_potential_dot_dot_acoustic,b_potential_acoustic,iphase)


! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,CPML_X_ONLY,CPML_Z_ONLY,IRIGHT,ILEFT,IBOTTOM,ITOP, &
    ZERO,ONE,TWO,TWO_THIRDS

  use specfem_par, only: nglob, &
                         ibool,ispec_is_acoustic, &
                         xix,xiz,gammax,gammaz,jacobian, &
                         rhostore, &
                         hprime_xx,hprimewgll_xx, &
                         hprime_zz,hprimewgll_zz,wxgll,wzgll, &
                         AXISYM,coord, is_on_the_axis,hprimeBar_xx,hprimeBarwglj_xx,xiglj,wxglj

  ! overlapping communication
  use specfem_par, only: nspec_inner_acoustic,nspec_outer_acoustic,phase_ispec_inner_acoustic

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob) :: b_potential_dot_dot_acoustic, b_potential_acoustic

  integer,intent(in) :: iphase

  ! local parameters
  integer :: ispec,i,j,k,iglob

  ! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,dux_dxl,dux_dzl
  real(kind=CUSTOM_REAL) :: xxi

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  ! material properties of the acoustic medium
  real(kind=CUSTOM_REAL) :: rhol

  integer :: num_elements,ispec_p

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

    ! acoustic spectral element
    if (.not. ispec_is_acoustic(ispec)) cycle

    ! first double loop over GLL points to compute and store gradients
    do j = 1,NGLLZ
      do i = 1,NGLLX
        ! derivative along x and along z
        dux_dxi = 0._CUSTOM_REAL
        dux_dgamma = 0._CUSTOM_REAL

        ! we can merge the two loops because NGLLX == NGLLZ
        if (AXISYM) then
          ! axisymmetric case
          if (is_on_the_axis(ispec)) then
            do k = 1,NGLLX
              dux_dxi = dux_dxi + b_potential_acoustic(ibool(k,j,ispec)) * hprimeBar_xx(i,k)
              dux_dgamma = dux_dgamma + b_potential_acoustic(ibool(i,k,ispec)) * hprime_zz(j,k)
            enddo
          else
            do k = 1,NGLLX
              dux_dxi = dux_dxi + b_potential_acoustic(ibool(k,j,ispec)) * hprime_xx(i,k)
              dux_dgamma = dux_dgamma + b_potential_acoustic(ibool(i,k,ispec)) * hprime_zz(j,k)
            enddo
          endif
        else
          ! default, non-axisymmetric case
          do k = 1,NGLLX
            dux_dxi = dux_dxi + b_potential_acoustic(ibool(k,j,ispec)) * hprime_xx(i,k)
            dux_dgamma = dux_dgamma + b_potential_acoustic(ibool(i,k,ispec)) * hprime_zz(j,k)
          enddo
        endif

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        ! derivatives of potential
        dux_dxl = dux_dxi * xixl + dux_dgamma * gammaxl
        dux_dzl = dux_dxi * xizl + dux_dgamma * gammazl

        if (AXISYM .and. is_on_the_axis(ispec) .and. i == 1) then ! dchi/dr=rho * u_r=0 on the axis
          dux_dxl = ZERO
        endif

        ! gets density
        rhol = rhostore(i,j,ispec)
        jacobianl = jacobian(i,j,ispec)

        if (AXISYM) then
          if (is_on_the_axis(ispec) .and. i == 1) then
            xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
            r_xiplus1(i,j) = xxi
          else if (is_on_the_axis(ispec)) then
            r_xiplus1(i,j) = coord(1,ibool(i,j,ispec))/(xiglj(i)+ONE)
          endif
        endif

        ! for acoustic medium also add integration weights
        if (AXISYM) then
          if (is_on_the_axis(ispec)) then
            tempx1(i,j) = wzgll(j) * r_xiplus1(i,j) * jacobianl * (xixl * dux_dxl + xizl * dux_dzl) / rhol
            tempx2(i,j) = wxglj(i) * r_xiplus1(i,j) * jacobianl * (gammaxl * dux_dxl + gammazl * dux_dzl) / rhol
          else
            tempx1(i,j) = wzgll(j) * coord(1,ibool(i,j,ispec)) * jacobianl * (xixl * dux_dxl + xizl * dux_dzl) / rhol
            tempx2(i,j) = wxgll(i) * coord(1,ibool(i,j,ispec)) * jacobianl * (gammaxl * dux_dxl + gammazl * dux_dzl) / rhol
          endif
        else
          tempx1(i,j) = wzgll(j) * jacobianl * (xixl * dux_dxl + xizl * dux_dzl) / rhol
          tempx2(i,j) = wxgll(i) * jacobianl * (gammaxl * dux_dxl + gammazl * dux_dzl) / rhol
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
        if (AXISYM) then
          ! axisymmetric case
          if (is_on_the_axis(ispec)) then
            do k = 1,NGLLX
              b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - &
                      (tempx1(k,j) * hprimeBarwglj_xx(k,i) + tempx2(i,k) * hprimewgll_zz(k,j))
            enddo
          else
            do k = 1,NGLLX
              b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - &
                      (tempx1(k,j) * hprimewgll_xx(k,i) + tempx2(i,k) * hprimewgll_zz(k,j))
            enddo
          endif
        else
          ! default, non-axisymmetric case
          do k = 1,NGLLX
            b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - &
                      (tempx1(k,j) * hprimewgll_xx(k,i) + tempx2(i,k) * hprimewgll_zz(k,j))
          enddo
        endif
      enddo ! second loop over the GLL points
    enddo

  enddo ! end of loop over all spectral elements

  end subroutine compute_forces_viscoacoustic_backward
