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


  subroutine compute_forces_acoustic_backward(b_minus_pressure_acoustic,b_minus_int_int_pressure_acoustic)


! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,CPML_X_ONLY,CPML_Z_ONLY,IRIGHT,ILEFT,IBOTTOM,ITOP, &
    ZERO,ONE,TWO,IEDGE1,IEDGE2,IEDGE3,IEDGE4

  use specfem_par, only: nglob,nspec, &
                         assign_external_model,ibool,kmato,ispec_is_acoustic, &
                         density,poroelastcoef,xix,xiz,gammax,gammaz,jacobian, &
                         vpext,rhoext, &
                         hprime_xx,hprimewgll_xx, &
                         hprime_zz,hprimewgll_zz,wxgll,wzgll, &
                         AXISYM,coord, is_on_the_axis,hprimeBar_xx,hprimeBarwglj_xx,xiglj,wxglj

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob) :: b_minus_pressure_acoustic, b_minus_int_int_pressure_acoustic

  ! local parameters
  integer :: ispec,i,j,k,iglob
  integer :: ifirstelem,ilastelem

  ! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,dux_dxl,dux_dzl
  real(kind=CUSTOM_REAL) :: xxi

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  ! material properties of the acoustic medium
  real(kind=CUSTOM_REAL) :: mul_relaxed,lambdal_relaxed,kappal,cpl,rhol

  ifirstelem = 1
  ilastelem = nspec


! loop over spectral elements
  do ispec = ifirstelem,ilastelem

    ! acoustic spectral element
    if (ispec_is_acoustic(ispec)) then

      rhol = density(1,kmato(ispec))

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX
          ! derivative along x and along z
          dux_dxi = 0._CUSTOM_REAL; dux_dgamma = 0._CUSTOM_REAL

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            if (AXISYM) then
              if (is_on_the_axis(ispec)) then
                dux_dxi = dux_dxi + b_minus_int_int_pressure_acoustic(ibool(k,j,ispec)) * hprimeBar_xx(i,k)
              else
                dux_dxi = dux_dxi + b_minus_int_int_pressure_acoustic(ibool(k,j,ispec)) * hprime_xx(i,k)
              endif
            else
              dux_dxi = dux_dxi + b_minus_int_int_pressure_acoustic(ibool(k,j,ispec)) * hprime_xx(i,k)
            endif
            dux_dgamma = dux_dgamma + b_minus_int_int_pressure_acoustic(ibool(i,k,ispec)) * hprime_zz(j,k)
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

          ! derivatives of the scalar field
          dux_dxl = dux_dxi * xixl + dux_dgamma * gammaxl
          dux_dzl = dux_dxi * xizl + dux_dgamma * gammazl

          if (AXISYM .and. is_on_the_axis(ispec) .and. i == 1) then ! dchi/dr=rho * u_r=0 on the axis
            dux_dxl = ZERO
          endif

          jacobianl = jacobian(i,j,ispec)

          ! if external density model
          if (assign_external_model) then
            rhol = rhoext(i,j,ispec)
          endif

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

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          if (assign_external_model) then
            rhol = rhoext(i,j,ispec)
            cpl = vpext(i,j,ispec)
            !assuming that in fluid(acoustic) part input cpl is defined by sqrt(kappal/rhol), &
            !which is not the same as in cpl input in elastic part
            kappal = rhol * cpl * cpl
          else
            lambdal_relaxed = poroelastcoef(1,1,kmato(ispec))
            mul_relaxed = poroelastcoef(2,1,kmato(ispec))
            kappal  = lambdal_relaxed + TWO * mul_relaxed/3._CUSTOM_REAL
            rhol = density(1,kmato(ispec))
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
            if (is_on_the_axis(ispec)) then
              do k = 1,NGLLX
                b_minus_pressure_acoustic(iglob) = b_minus_pressure_acoustic(iglob) - &
                        (tempx1(k,j) * hprimeBarwglj_xx(k,i) + tempx2(i,k) * hprimewgll_zz(k,j))
              enddo
            else
              do k = 1,NGLLX
                b_minus_pressure_acoustic(iglob) = b_minus_pressure_acoustic(iglob) - &
                        (tempx1(k,j) * hprimewgll_xx(k,i) + tempx2(i,k) * hprimewgll_zz(k,j))
              enddo
            endif
          else
            do k = 1,NGLLX
              b_minus_pressure_acoustic(iglob) = b_minus_pressure_acoustic(iglob) - &
                        (tempx1(k,j) * hprimewgll_xx(k,i) + tempx2(i,k) * hprimewgll_zz(k,j))
            enddo
          endif
        enddo ! second loop over the GLL points
      enddo

    endif ! end of test if acoustic element
  enddo ! end of loop over all spectral elements

  end subroutine compute_forces_acoustic_backward
