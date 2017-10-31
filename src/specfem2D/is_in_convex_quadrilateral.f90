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

  subroutine is_in_convex_quadrilateral(elmnt_coords, x_coord, z_coord, is_in)

  implicit none

  double precision, dimension(2,4)  :: elmnt_coords
  double precision, intent(in)  :: x_coord, z_coord
  logical, intent(out)  :: is_in

  real :: x1, x2, x3, x4, z1, z2, z3, z4
  real  :: normal1, normal2, normal3, normal4

  x1 = elmnt_coords(1,1)
  x2 = elmnt_coords(1,2)
  x3 = elmnt_coords(1,3)
  x4 = elmnt_coords(1,4)
  z1 = elmnt_coords(2,1)
  z2 = elmnt_coords(2,2)
  z3 = elmnt_coords(2,3)
  z4 = elmnt_coords(2,4)

  normal1 = (z_coord-z1) * (x2-x1) - (x_coord-x1) * (z2-z1)
  normal2 = (z_coord-z2) * (x3-x2) - (x_coord-x2) * (z3-z2)
  normal3 = (z_coord-z3) * (x4-x3) - (x_coord-x3) * (z4-z3)
  normal4 = (z_coord-z4) * (x1-x4) - (x_coord-x4) * (z1-z4)

  if ((normal1 < 0) .or. (normal2 < 0) .or. (normal3 < 0) .or. (normal4 < 0)) then
    is_in = .false.
  else
    is_in = .true.
  endif

  end subroutine is_in_convex_quadrilateral

