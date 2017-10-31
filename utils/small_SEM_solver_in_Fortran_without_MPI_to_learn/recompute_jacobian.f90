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

! Recompute 2D jacobian at a given point in a 4-node or 9-node element
! Compute also the global coordinates of the point defined by: (xi,gamma,ispec)

  subroutine recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian,coorg,knods,ispec,ngnod,nspec,npgeo,NDIM)

  implicit none

  integer, intent(in) :: ispec,ngnod,nspec,npgeo,NDIM
  double precision, intent(out) :: x,z,xix,xiz,gammax,gammaz,jacobian
  double precision, intent(in) :: xi,gamma

  integer, intent(in) :: knods(ngnod,nspec)
  double precision, intent(in) :: coorg(NDIM,npgeo)

! 2D shape functions and their derivatives at receiver
  double precision shape2D(ngnod)
  double precision dershape2D(NDIM,ngnod)

  double precision xxi,zxi,xgamma,zgamma,xelm,zelm

  integer ia,nnum

! recompute jacobian for any (xi,gamma) point, not necessarily a GLL point

! create the 2D shape functions and then the Jacobian
  call define_shape_functions(shape2D,dershape2D,xi,gamma,ngnod,NDIM)

! compute coordinates and jacobian matrix
  x = 0.d0
  z = 0.d0

  xxi = 0.d0
  zxi = 0.d0
  xgamma = 0.d0
  zgamma = 0.d0

  do ia=1,ngnod

    nnum = knods(ia,ispec)

    xelm = coorg(1,nnum)
    zelm = coorg(2,nnum)

    x = x + shape2D(ia)*xelm
    z = z + shape2D(ia)*zelm

    xxi = xxi + dershape2D(1,ia)*xelm
    zxi = zxi + dershape2D(1,ia)*zelm
    xgamma = xgamma + dershape2D(2,ia)*xelm
    zgamma = zgamma + dershape2D(2,ia)*zelm

  enddo

  jacobian = xxi*zgamma - xgamma*zxi

! the Jacobian is negative, so far this means that there is an error in the mesh
! therefore print the coordinates of the mesh points of this element
  if (jacobian <= 0.d0) then
! print the coordinates of the mesh points of this element
    print *, 'ispec = ', ispec
    print *, 'ngnod = ', ngnod
    do ia=1,ngnod
      nnum = knods(ia,ispec)
      xelm = coorg(1,nnum)
      zelm = coorg(2,nnum)
      print *,'node ', ia,' x,y = ',xelm,zelm
    enddo
    stop 'error: negative 2D Jacobian found'
  endif

! invert the relation
  xix = zgamma / jacobian
  gammax = - zxi / jacobian
  xiz = - xgamma / jacobian
  gammaz = xxi / jacobian

  end subroutine recompute_jacobian

