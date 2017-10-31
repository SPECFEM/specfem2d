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

  subroutine enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic)

! free surface for an acoustic medium
! if acoustic, the free surface condition is a Dirichlet condition for the potential,
! not Neumann, in order to impose zero pressure at the surface

  use constants, only: CUSTOM_REAL,ZERO

  use specfem_par, only: acoustic_surface,ibool,nelem_acoustic_surface,nglob,this_ibool_is_a_periodic_edge,coord

  implicit none

!! DK DK added this to exclude the bottom surface of the mesh for Laurent Guillon in order not to create
!! DK DK a Dirichlet condition for the potential, i.e. in order to create a Neumann condition for the potential
!! DK DK i.e. in order to create a rigid bottom surface instead of a free surface; I test the Z coordinate of the mesh point
  logical, parameter :: ENFORCE_RIGID_SURFACE_BOTTOM = .false.
! this should be a bit bigger than the Z coordinate of the bottom; since it is Z = 0 in Laurent Guillon's test
! we use this small value (its actual value does not matter as long as it is smaller than the Z coordinate of the top surface)
  double precision, parameter :: Zlimit = 0.0000001d0

  real(kind=CUSTOM_REAL), dimension(nglob) :: potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic

!---
!--- local variables
!---

  integer :: ispec_acoustic_surface,ispec,i,j,iglob

  ! checks if anything to do
  if (nelem_acoustic_surface == 0) return

  do ispec_acoustic_surface = 1, nelem_acoustic_surface

    ispec = acoustic_surface(1,ispec_acoustic_surface)

    do j = acoustic_surface(4,ispec_acoustic_surface), acoustic_surface(5,ispec_acoustic_surface)
      do i = acoustic_surface(2,ispec_acoustic_surface), acoustic_surface(3,ispec_acoustic_surface)
        iglob = ibool(i,j,ispec)
        ! make sure that an acoustic free surface is not enforced on periodic edges
        if (.not. this_ibool_is_a_periodic_edge(iglob)) then
!! DK DK added this to exclude the bottom surface of the mesh for Laurent Guillon in order not to create
!! DK DK a Dirichlet condition for the potential, i.e. in order to create a Neumann condition for the potential
!! DK DK i.e. in order to create a rigid bottom surface instead of a free surface; I test the Z coordinate of the mesh point
        if (.not. ENFORCE_RIGID_SURFACE_BOTTOM .or. coord(2,iglob) >= Zlimit) then
          potential_acoustic(iglob) = ZERO
          potential_dot_acoustic(iglob) = ZERO
          potential_dot_dot_acoustic(iglob) = ZERO
        endif
        endif
      enddo
    enddo

  enddo

  end subroutine enforce_acoustic_free_surface

