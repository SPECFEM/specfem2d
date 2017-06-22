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

  subroutine compute_curl_one_element(ispec,curl_element)

! compute curl in (poro)elastic elements (for rotational study)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,ZERO

  use specfem_par, only: displ_elastic,displs_poroelastic,ispec_is_elastic,ispec_is_poroelastic, &
                         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,myrank

  implicit none

  integer,intent(in) :: ispec

  ! curl in an element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(out) :: curl_element

  ! local parameters
  ! jacobian
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl

  ! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: duz_dxl,dux_dzl
  integer :: i,j,k

  ! initializes
  curl_element(:,:) = 0._CUSTOM_REAL

  if (ispec_is_elastic(ispec)) then
    ! elastic element
    do j = 1,NGLLZ
      do i = 1,NGLLX
        ! derivative along x and along z
        dux_dxi = ZERO
        duz_dxi = ZERO

        dux_dgamma = ZERO
        duz_dgamma = ZERO

        ! first double loop over GLL points to compute and store gradients
        ! we can merge the two loops because NGLLX == NGLLZ
        do k = 1,NGLLX
          dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
          duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
          dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
          duz_dgamma = duz_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        ! derivatives of displacement
        dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl
        duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl

        ! store pressure
        curl_element(i,j) = - 0.5d0 * (dux_dzl - duz_dxl)
      enddo
    enddo

  else if (ispec_is_poroelastic(ispec)) then
    ! poro-elastic element
    do j = 1,NGLLZ
      do i = 1,NGLLX
        ! derivative along x and along z
        dux_dxi = ZERO
        duz_dxi = ZERO

        dux_dgamma = ZERO
        duz_dgamma = ZERO

        ! first double loop over GLL points to compute and store gradients
        ! we can merge the two loops because NGLLX == NGLLZ
        do k = 1,NGLLX
          dux_dxi = dux_dxi + displs_poroelastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
          duz_dxi = duz_dxi + displs_poroelastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
          dux_dgamma = dux_dgamma + displs_poroelastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
          duz_dgamma = duz_dgamma + displs_poroelastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        ! derivatives of displacement
        dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl
        duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl

        ! store pressure
        curl_element(i,j) = - 0.5d0 * (dux_dzl - duz_dxl)
      enddo
    enddo

  else
    ! safety stop
    call exit_MPI(myrank,'no curl in acoustic')

  endif ! end of test if acoustic or elastic element

  end subroutine compute_curl_one_element

