
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
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and Inria at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

! Recompute 2D jacobian at a given point in a 4-node or 9-node element
! Compute also the global coordinates of the point defined by: (xi,gamma,ispec)

  subroutine recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian,coorg,knods,ispec,ngnod,nspec,npgeo, &
                      stop_if_negative_jacobian)

  implicit none

  include "constants.h"

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
! and also create an OpenDX file to visualize it
  if(jacobian <= ZERO .and. stop_if_negative_jacobian) then

! print the coordinates of the mesh points of this element
    print *, 'ispec = ', ispec
    print *, 'ngnod = ', ngnod
    do ia=1,ngnod
      nnum = knods(ia,ispec)
      xelm = coorg(1,nnum)
      zelm = coorg(2,nnum)
      print *,'node ', ia,' x,y = ',xelm,zelm
    enddo

! create an OpenDX file to visualize this element
    open(unit=11,file='DX_first_element_with_negative_jacobian.dx',status='unknown')

! output the points (the mesh is flat therefore the third coordinate is zero)
    write(11,*) 'object 1 class array type float rank 1 shape 3 items ',ngnod,' data follows'
    do ia=1,ngnod
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

    call exit_MPI('negative 2D Jacobian, element saved in DX_first_element_with_negative_jacobian.dx')
  endif

! invert the relation
  xix = zgamma / jacobian
  gammax = - zxi / jacobian
  xiz = - xgamma / jacobian
  gammaz = xxi / jacobian

  end subroutine recompute_jacobian

