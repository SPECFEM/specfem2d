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

! From array 'surface' (element, type : node/edge, node(s) ) that describes the
! elastic fixed surface, determines the points (ixmin, ixmax, izmin and izmax) on the surface for each element.
! We chose to have ixmin <= ixmax and izmin <= izmax, so as to be able to have DO loops on it with an increment of +1.

  subroutine construct_elastic_fixed_surface ()

  use specfem_par, only: ngnod, knods, nelem_elastic_fixed_surface, elastic_fixed_edges, elastic_fixed_surface

  implicit none

  integer  :: i, k
  integer  :: ixmin, ixmax
  integer  :: izmin, izmax
  integer, dimension(ngnod)  :: n
  integer  :: e1, e2
  integer  :: type_elastic

  do i = 1, nelem_elastic_fixed_surface
     elastic_fixed_surface(1,i) = elastic_fixed_edges(1,i) ! Here we do a copy
     type_elastic = elastic_fixed_edges(2,i)
     e1 = elastic_fixed_edges(3,i)
     e2 = elastic_fixed_edges(4,i)
     do k = 1, ngnod
        n(k) = knods(k,elastic_fixed_surface(1,i))
     enddo

     call get_elastic_fixed_edge ( ngnod, n, type_elastic, e1, e2, ixmin, ixmax, izmin, izmax )

     elastic_fixed_surface(2,i) = ixmin
     elastic_fixed_surface(3,i) = ixmax
     elastic_fixed_surface(4,i) = izmin
     elastic_fixed_surface(5,i) = izmax

  enddo

  end subroutine construct_elastic_fixed_surface

!-----------------------------------------------
! Get the points (ixmin, ixmax, izmin and izmax) on an node/edge for one element.
!-----------------------------------------------
  subroutine get_elastic_fixed_edge ( ngnod, n, type, e1, e2, ixmin, ixmax, izmin, izmax )

  use constants, only: NGLLX,NGLLZ

  implicit none

  integer, intent(in)  :: ngnod
  integer, dimension(ngnod), intent(in)  :: n
  integer, intent(in)  :: type, e1, e2
  integer, intent(out)  :: ixmin, ixmax, izmin, izmax


  if (type == 1) then
     if (e1 == n(1)) then
        ixmin = 1
        ixmax = 1
        izmin = 1
        izmax = 1
     endif
     if (e1 == n(2)) then
        ixmin = NGLLX
        ixmax = NGLLX
        izmin = 1
        izmax = 1
     endif
     if (e1 == n(3)) then
        ixmin = NGLLX
        ixmax = NGLLX
        izmin = NGLLZ
        izmax = NGLLZ
     endif
     if (e1 == n(4)) then
        ixmin = 1
        ixmax = 1
        izmin = NGLLZ
        izmax = NGLLZ
     endif

  else
     if (e1 == n(1)) then
        ixmin = 1
        izmin = 1
        if (e2 == n(2)) then
           ixmax = NGLLX
           izmax = 1

        endif
        if (e2 == n(4)) then
           ixmax = 1
           izmax = NGLLZ

        endif
     endif
     if (e1 == n(2)) then
        ixmin = NGLLX
        izmin = 1
        if (e2 == n(3)) then
           ixmax = NGLLX
           izmax = NGLLZ

        endif
        if (e2 == n(1)) then
           ixmax = ixmin
           ixmin = 1
           izmax = 1

        endif
     endif
     if (e1 == n(3)) then
        ixmin = NGLLX
        izmin = NGLLZ
        if (e2 == n(4)) then
           ixmax = ixmin
           ixmin = 1
           izmax = NGLLZ

        endif
        if (e2 == n(2)) then
           ixmax = NGLLX
           izmax = izmin
           izmin = 1

        endif
     endif
     if (e1 == n(4)) then
        ixmin = 1
        izmin = NGLLZ
        if (e2 == n(1)) then
           ixmax = 1
           izmax = izmin
           izmin = 1

        endif
        if (e2 == n(3)) then
           ixmax = NGLLX
           izmax = NGLLZ

        endif
     endif
  endif

  end subroutine get_elastic_fixed_edge


!-----------------------------------------------
! set the accel_elastic, veloc_elastic, and displ_elastic on the points (ixmin, ixmax, izmin and izmax)
! that on an node/edge for one element to be zero
!-----------------------------------------------

  subroutine enforce_elastic_fixed_surface(accel_elastic,veloc_elastic,displ_elastic)

  use constants, only: CUSTOM_REAL,ZERO,NDIM

  use specfem_par, only: elastic_fixed_surface,ibool,nelem_elastic_fixed_surface,nglob, &
                         this_ibool_is_a_periodic_edge

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob) :: accel_elastic,veloc_elastic,displ_elastic

!---
!--- local variables
!---

  integer :: ispec_elastic_fixed_surface,ispec,i,j,iglob

  do ispec_elastic_fixed_surface = 1, nelem_elastic_fixed_surface

    ispec = elastic_fixed_surface(1,ispec_elastic_fixed_surface)

    do j = elastic_fixed_surface(4,ispec_elastic_fixed_surface), elastic_fixed_surface(5,ispec_elastic_fixed_surface)
      do i = elastic_fixed_surface(2,ispec_elastic_fixed_surface), elastic_fixed_surface(3,ispec_elastic_fixed_surface)
        iglob = ibool(i,j,ispec)
        ! make sure that an acoustic free surface is not enforced on periodic edges
        if (.not. this_ibool_is_a_periodic_edge(iglob)) then
          displ_elastic(:,iglob) = ZERO
          veloc_elastic(:,iglob) = ZERO
          accel_elastic(:,iglob) = ZERO
        endif
      enddo
    enddo

  enddo

  end subroutine enforce_elastic_fixed_surface
