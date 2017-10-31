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


! *******************
! meshing subroutines
! *******************

!--- global node number

  integer function num(i,j,nx)

  implicit none

  integer i,j,nx

  num = j*(nx+1) + i + 1

  end function num


!---  global node number (when ngnod==4).
  integer function num_4(i,j,nx)

  implicit none

  integer i,j,nx

  num_4 = j*(nx+1) + i + 1

  end function num_4


!---  global node number (when ngnod==9).
  integer function num_9(i,j,nx,nz)

  implicit none

  integer i,j,nx,nz


  if ((mod(i,2) == 0) .and. (mod(j,2) == 0)) then
     num_9 = j/2 * (nx+1) + i/2 + 1
  else
     if (mod(j,2) == 0) then
        num_9 = (nx+1)*(nz+1) + j/2 * nx + ceiling(real(i)/real(2))
     else
        num_9 = (nx+1)*(nz+1) + nx*(nz+1) + floor(real(j)/real(2))*(nx*2+1) + i + 1

     endif
  endif

  end function num_9
