
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
!                         Dimitri Komatitsch
!                     University of Pau, France
!
!                          (c) April 2007
!
!========================================================================

! recompute 2D jacobian at a given point in a 4-node or 9-node element

  subroutine recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian,coorg,knods,ispec,ngnod,nspec,npgeo)

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

! recompute jacobian for any (xi,gamma) point, not necessarily a GLL point

! create the 2D shape functions and the Jacobian
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

  if(jacobian <= ZERO) stop '2D Jacobian undefined'

! invert the relation
  xix = zgamma / jacobian
  gammax = - zxi / jacobian
  xiz = - xgamma / jacobian
  gammaz = xxi / jacobian

  end subroutine recompute_jacobian

