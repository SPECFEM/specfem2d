
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.0
!                   ------------------------------
!
!                         Dimitri Komatitsch
!          Universite de Pau et des Pays de l'Adour, France
!
!                          (c) May 2004
!
!========================================================================

  subroutine q49spec(shapeint,dershape,xix,xiz,gammax,gammaz,jacobian,xi, &
              coorg,knods,ngnod,nspec,npgeo, &
              xirec,etarec,flagrange,iptsdisp)

!=======================================================================
!
!                       set up the jacobian matrix
!                       for the isoparametric transformation of the
!                       spectral macroblocs.
!                       The routine can handle 4 or 9 control nodes.
!                       The control nodes are defined as follows:
!
!                               4 . . . . 7 . . . . 3
!                               .                   .
!                               .         t         .
!                               .                   .
!                               8         9  s      6
!                               .                   .
!                               .                   .
!                               .                   .
!                               1 . . . . 5 . . . . 2
!
!                           Local coordinate system : s,t
!
!=======================================================================

  implicit none

  include "constants.h"

  integer ngnod,nspec,npgeo,iptsdisp

  integer knods(ngnod,nspec)
  double precision shapeint(ngnod,iptsdisp,iptsdisp)
  double precision dershape(NDIME,ngnod,NGLLX,NGLLX)
  double precision coorg(NDIME,npgeo)
  double precision xi(NGLLX)
  double precision xirec(iptsdisp),etarec(iptsdisp)
  double precision flagrange(NGLLX,iptsdisp)

  double precision xix(NGLLX,NGLLZ,nspec)
  double precision xiz(NGLLX,NGLLZ,nspec)
  double precision gammax(NGLLX,NGLLZ,nspec)
  double precision gammaz(NGLLX,NGLLZ,nspec)
  double precision jacobian(NGLLX,NGLLZ,nspec)

  integer l1,l2,ispec,in,nnum,ip1,ip2,i,k
  double precision s,sp,sm,t,tp,tm,s2,t2,ss,tt,st
  double precision xjac2_11,xjac2_21,xjac2_12,xjac2_22

  double precision, external :: hgll

!
!----    compute the jacobian matrix at the integration points
!

  do ispec = 1,nspec

  do ip1 = 1,NGLLX
  do ip2 = 1,NGLLZ

    xjac2_11 = ZERO
    xjac2_21 = ZERO
    xjac2_12 = ZERO
    xjac2_22 = ZERO

    do in = 1,ngnod

      nnum = knods(in,ispec)

      xjac2_11 = xjac2_11 + dershape(1,in,ip1,ip2)*coorg(1,nnum)
      xjac2_21 = xjac2_21 + dershape(1,in,ip1,ip2)*coorg(2,nnum)
      xjac2_12 = xjac2_12 + dershape(2,in,ip1,ip2)*coorg(1,nnum)
      xjac2_22 = xjac2_22 + dershape(2,in,ip1,ip2)*coorg(2,nnum)

    enddo

!
!----    invert the jacobian matrix
!

    jacobian(ip1,ip2,ispec) = xjac2_11*xjac2_22 - xjac2_12*xjac2_21

    if(jacobian(ip1,ip2,ispec) <= zero) stop 'Error: Jacobian undefined!'

    xix(ip1,ip2,ispec) = xjac2_22 / jacobian(ip1,ip2,ispec)
    gammax(ip1,ip2,ispec) = - xjac2_21 / jacobian(ip1,ip2,ispec)
    xiz(ip1,ip2,ispec) = - xjac2_12 / jacobian(ip1,ip2,ispec)
    gammaz(ip1,ip2,ispec) = xjac2_11 / jacobian(ip1,ip2,ispec)

  enddo
  enddo

  enddo

!---- calcul des coordonnees interpolees avec les fonctions de forme
!---- interpolation sur grille reguliere en (xi,eta)

  do i=1,iptsdisp
    xirec(i)  = 2.d0*dble(i-1)/dble(iptsdisp-1) - 1.d0
    etarec(i) = 2.d0*dble(i-1)/dble(iptsdisp-1) - 1.d0
  enddo

!---- calcul des interpolateurs de Lagrange (suppose NGLLX = NGLLZ)
  do i=1,NGLLX
    do k=1,iptsdisp
      flagrange(i,k) = hgll(i-1,xirec(k),xi,NGLLX)
    enddo
  enddo

!
!---- set up the shape functions for the interpolated grid
!
  if(ngnod == 4) then
!
!----    4-noded rectangular element
!
 do l2 = 1,iptsdisp

    t  = etarec(l2)

    do l1 = 1,iptsdisp

       s  = xirec(l1)

       sp = s + one
       sm = s - one
       tp = t + one
       tm = t - one

!
!----          corner nodes
!
       shapeint(1,l1,l2) = quart * sm * tm
       shapeint(2,l1,l2) = - quart * sp * tm
       shapeint(3,l1,l2) = quart * sp * tp
       shapeint(4,l1,l2) = - quart * sm * tp

    enddo
 enddo

  else if(ngnod == 9) then
!
!----    9-noded rectangular element
!
 do l2 = 1,iptsdisp

    t  = etarec(l2)

    do l1 = 1,iptsdisp

       s  = xirec(l1)

       sp = s + one
       sm = s - one
       tp = t + one
       tm = t - one
       s2 = s * two
       t2 = t * two
       ss = s * s
       tt = t * t
       st = s * t

!
!----          corner nodes
!
       shapeint(1,l1,l2) = quart * sm * st * tm
       shapeint(2,l1,l2) = quart * sp * st * tm
       shapeint(3,l1,l2) = quart * sp * st * tp
       shapeint(4,l1,l2) = quart * sm * st * tp

!
!----          midside nodes
!
       shapeint(5,l1,l2) = half * tm * t * (one - ss)
       shapeint(6,l1,l2) = half * sp * s * (one - tt)
       shapeint(7,l1,l2) = half * tp * t * (one - ss)
       shapeint(8,l1,l2) = half * sm * s * (one - tt)

!
!----          center node
!
       shapeint(9,l1,l2) = (one - ss) * (one - tt)

    enddo
 enddo

  endif

  end subroutine q49spec

