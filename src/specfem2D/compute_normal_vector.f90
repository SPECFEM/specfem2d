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

  subroutine compute_normal_vector( angle, n1_x, n2_x, n3_x, n4_x, n1_z, n2_z, n3_z, n4_z )

  use constants, only: TINYVAL,PI

  implicit none

  double precision :: angle
  double precision :: n1_x, n2_x, n3_x, n4_x, n1_z, n2_z, n3_z, n4_z

  double precision  :: theta1, theta2, theta3
  double precision  :: costheta1, costheta2, costheta3

  if (abs(n2_z - n1_z) < TINYVAL) then
     costheta1 = 0
  else
     costheta1 = (n2_z - n1_z) / sqrt((n2_x - n1_x)**2 + (n2_z - n1_z)**2)
  endif
  if (abs(n3_z - n2_z) < TINYVAL) then
     costheta2 = 0
  else
     costheta2 = (n3_z - n2_z) / sqrt((n3_x - n2_x)**2 + (n3_z - n2_z)**2)
  endif
  if (abs(n4_z - n3_z) < TINYVAL) then
     costheta3 = 0
  else
    costheta3 = (n4_z - n3_z) / sqrt((n4_x - n3_x)**2 + (n4_z - n3_z)**2)
  endif

  theta1 = - sign(1.d0,n2_x - n1_x) * acos(costheta1)
  theta2 = - sign(1.d0,n3_x - n2_x) * acos(costheta2)
  theta3 = - sign(1.d0,n4_x - n3_x) * acos(costheta3)

  ! a sum is needed here because in the case of a source force vector
  ! users can give an angle with respect to the normal to the topography surface,
  ! in which case we must compute the normal to the topography
  ! and add it the existing rotation angle
  angle = angle + (theta1 + theta2 + theta3) / 3.d0 + PI/2.d0

  end subroutine compute_normal_vector


!
!-------------------------------------------------------------------------------------------------
!

  subroutine tri_quad(n, n1, nnodes)

  implicit none

  integer  :: n1, nnodes
  integer, dimension(4)  :: n


  n(2) = n1

  if (n1 == 1) then
     n(1) = nnodes
  else
     n(1) = n1-1
  endif

  if (n1 == nnodes) then
     n(3) = 1
  else
     n(3) = n1+1
  endif

  if (n(3) == nnodes) then
     n(4) = 1
  else
     n(4) = n(3)+1
  endif


  end subroutine tri_quad

