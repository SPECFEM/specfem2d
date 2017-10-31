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

! recompute 2D jacobian at a given point in a 4-node or 9-node element

  subroutine recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian,coord_of_anchor_points,ispec,nspec)

  implicit none

  include "constants.h"

  integer ispec,nspec
  double precision x,z,xix,xiz,gammax,gammaz
  double precision xi,gamma,jacobian

! geometrical coordinates of the anchor points of the elements of the mesh
  double precision, dimension(NDIM,NGNOD,NSPEC) :: coord_of_anchor_points

! 2D shape functions and their derivatives at receiver
  double precision shape2D(ngnod)
  double precision dershape2D(NDIM,ngnod)

  double precision xxi,zxi,xgamma,zgamma,xelm,zelm

  integer ia

! recompute jacobian for any (xi,gamma) point, not necessarily a GLL point

! create the 2D shape functions and then the Jacobian
  call define_shape_functions(shape2D,dershape2D,xi,gamma)

! compute coordinates and jacobian matrix
  x = ZERO
  z = ZERO

  xxi = ZERO
  zxi = ZERO
  xgamma = ZERO
  zgamma = ZERO

  do ia=1,ngnod
    xelm = coord_of_anchor_points(1,ia,ispec)
    zelm = coord_of_anchor_points(2,ia,ispec)

    x = x + shape2D(ia)*xelm
    z = z + shape2D(ia)*zelm

    xxi = xxi + dershape2D(1,ia)*xelm
    zxi = zxi + dershape2D(1,ia)*zelm
    xgamma = xgamma + dershape2D(2,ia)*xelm
    zgamma = zgamma + dershape2D(2,ia)*zelm
  enddo

  jacobian = xxi*zgamma - xgamma*zxi

! if the Jacobian is negative there is an error in the mesh
  if (jacobian <= ZERO) stop 'error: negative Jacobian found, the mesh is not correct'

! invert the relation
  xix = zgamma / jacobian
  gammax = - zxi / jacobian
  xiz = - xgamma / jacobian
  gammaz = xxi / jacobian

  end subroutine recompute_jacobian

