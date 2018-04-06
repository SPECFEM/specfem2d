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

! Recompute 2D jacobian at a given point in a 4-node or 9-node element.
! Compute also the global coordinates of the point defined by: (xi,gamma,ispec).

  subroutine recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian,coorg,knods,ispec,ngnod,nspec,npgeo, &
                      stop_if_negative_jacobian)

  use constants, only: NDIM,ZERO

  implicit none

  integer ispec,ngnod,nspec,npgeo
  double precision x,z,xix,xiz,gammax,gammaz
  double precision xi,gamma,jacobian

  integer knods(ngnod,nspec)
  double precision coorg(NDIM,npgeo)

! 2D shape functions and their derivatives at receiver
  double precision shape2D(ngnod)
  double precision dershape2D(NDIM,ngnod)

  double precision xxi,zxi,xgamma,zgamma,xelm,zelm

  integer ia,nnum

  logical stop_if_negative_jacobian

! only one problematic element is output to OpenDX for now in case of elements with a negative Jacobian
  integer, parameter :: ntotspecAVS_DX = 1

! recompute jacobian for any (xi,gamma) point, not necessarily a GLL point

! create the 2D shape functions and then the Jacobian
  call define_shape_functions(shape2D,dershape2D,xi,gamma,ngnod)

! compute coordinates and jacobian matrix
  x = ZERO
  z = ZERO

  xxi = ZERO
  zxi = ZERO
  xgamma = ZERO
  zgamma = ZERO

  do ia = 1,ngnod

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
! and also create an OpenDX file to visualize it
  if (jacobian <= ZERO .and. stop_if_negative_jacobian) then

! print the coordinates of the mesh points of this element
    print *, 'ispec = ', ispec
    print *, 'ngnod = ', ngnod
    do ia = 1,ngnod
      nnum = knods(ia,ispec)
      xelm = coorg(1,nnum)
      zelm = coorg(2,nnum)
      print *,'node ', ia,' x,y = ',xelm,zelm
    enddo

! create an OpenDX file to visualize this element
    open(unit=11,file='DX_first_element_with_negative_jacobian.dx',status='unknown')

! output the points (the mesh is flat therefore the third coordinate is zero)
    write(11,*) 'object 1 class array type float rank 1 shape 3 items ',ngnod,' data follows'
    do ia = 1,ngnod
      nnum = knods(ia,ispec)
      xelm = coorg(1,nnum)
      zelm = coorg(2,nnum)
      write(11,*) xelm,zelm,' 0'
    enddo

! output the element (use its four corners only for now)
    write(11,*) 'object 2 class array type int rank 1 shape 4 items ',ntotspecAVS_DX,' data follows'
! point order in OpenDX is 1,4,2,3 *not* 1,2,3,4 as in AVS
    write(11,*) '0 3 1 2'

! output element data
    write(11,*) 'attribute "element type" string "quads"'
    write(11,*) 'attribute "ref" string "positions"'
    write(11,*) 'object 3 class array type float rank 0 items ',ntotspecAVS_DX,' data follows'

! output dummy data value
    write(11,*) '1'

! define OpenDX field
    write(11,*) 'attribute "dep" string "connections"'
    write(11,*) 'object "irregular positions irregular connections" class field'
    write(11,*) 'component "positions" value 1'
    write(11,*) 'component "connections" value 2'
    write(11,*) 'component "data" value 3'
    write(11,*) 'end'

! close OpenDX file
    close(11)

    call stop_the_code('negative 2D Jacobian, element saved in DX_first_element_with_negative_jacobian.dx')
  endif

! invert the relation
  xix = zgamma / jacobian
  gammax = - zxi / jacobian
  xiz = - xgamma / jacobian
  gammaz = xxi / jacobian

  end subroutine recompute_jacobian

