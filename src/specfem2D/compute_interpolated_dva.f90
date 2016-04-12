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
!=====================================================================

  subroutine compute_interpolated_dva(irec,ispec,vector_field_element,pressure_element,curl_element, &
                                      vx,vz,vcurl)

  use constants,only: CUSTOM_REAL,NDIM,NGLLX,NGLLZ
  use specfem_par,only: hxir_store,hgammar_store,ibool,seismotype

  implicit none

  integer,intent(in) :: irec,ispec
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ),intent(in) :: vector_field_element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(in) :: pressure_element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(in) :: curl_element

  double precision,intent(out) :: vx,vz,vcurl

  ! local parameters
  double precision :: dxd,dzd,dcurld,hlagrange
  integer :: i,j,iglob

  ! initializes interpolated values
  vx = 0.d0
  vz = 0.d0
  vcurl = 0.d0

  do j = 1,NGLLZ
    do i = 1,NGLLX
      iglob = ibool(i,j,ispec)
      hlagrange = hxir_store(irec,i)*hgammar_store(irec,j)

      ! displacement/velocity/acceleration/pressure value (depending on seismotype)
      select case (seismotype)
      case (1,2,3)
        ! displacement/velocity/acceleration
        dxd = vector_field_element(1,i,j)
        dzd = vector_field_element(2,i,j)
        ! computes interpolated field
        vx = vx + dxd*hlagrange
        vz = vz + dzd*hlagrange

      case (4,6)
        ! pressure
        dxd = pressure_element(i,j)
        ! computes interpolated field
        vx = vx + dxd*hlagrange

      case (5)
        ! displacement
        dxd = vector_field_element(1,i,j)
        dzd = vector_field_element(2,i,j)
        ! curl of displacement
        dcurld = curl_element(i,j)
        ! computes interpolated field
        vx = vx + dxd*hlagrange
        vz = vz + dzd*hlagrange
        vcurl = vcurl + dcurld*hlagrange

      case default
        stop 'Invalid seismotype for writing seismograms'
      end select
    enddo
  enddo

  end subroutine compute_interpolated_dva


