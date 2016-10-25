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

  subroutine enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic)

! free surface for an acoustic medium
! if acoustic, the free surface condition is a Dirichlet condition for the potential,
! not Neumann, in order to impose zero pressure at the surface

  use constants, only: CUSTOM_REAL,ZERO

  use specfem_par, only: acoustic_surface,ibool,nelem_acoustic_surface,nglob,this_ibool_is_a_periodic_edge

  implicit none

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
          potential_acoustic(iglob) = ZERO
          potential_dot_acoustic(iglob) = ZERO
          potential_dot_dot_acoustic(iglob) = ZERO
        endif
      enddo
    enddo

  enddo

  end subroutine enforce_acoustic_free_surface

