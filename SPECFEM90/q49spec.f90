
!=====================================================================
!
!                 S p e c f e m  V e r s i o n  4 . 2
!                 -----------------------------------
!
!                         Dimitri Komatitsch
!    Department of Earth and Planetary Sciences - Harvard University
!                         Jean-Pierre Vilotte
!                 Departement de Sismologie - IPGP - Paris
!                           (c) June 1998
!
!=====================================================================

  subroutine q49spec(shapeint,dershape,dvolu,xjaci,xi, &
              coorg,knods,ngnod,nxgll,nygll,ndime,nspec,npgeo, &
              xirec,etarec,flagrange,iptsdisp)
!
!=======================================================================
!
!     "q 4 9 s p e c" : set up the jacobian matrix
!                       for the isoparametric transformation of the
!                       spectral macroblocs.
!                       The routine is able to deal with
!                       4 or 9 control nodes for each bloc.
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
!

  implicit none

  integer ngnod,nxgll,nygll,ndime,nspec,npgeo,iptsdisp

  integer knods(ngnod,nspec)
  double precision shapeint(ngnod,iptsdisp,iptsdisp)
  double precision dershape(ndime,ngnod,nxgll,nxgll)
  double precision dvolu(nspec,nxgll,nxgll)
  double precision xjaci(nspec,ndime,ndime,nxgll,nxgll)
  double precision coorg(ndime,npgeo)
  double precision xi(nxgll)
  double precision xirec(iptsdisp),etarec(iptsdisp)
  double precision flagrange(0:nxgll-1,iptsdisp)

  double precision, parameter :: &
       zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,quart=0.25d0

  integer l1,l2,ispel,in,nnum,ip1,ip2,i,k
  double precision s,sp,sm,t,tp,tm,s2,t2,ss,tt,st
  double precision xjac2_11,xjac2_21,xjac2_12,xjac2_22

  double precision, external :: hgll

!
!-----------------------------------------------------------------------
!

!
!----    compute the jacobian matrix at the integration points
!

  do ispel = 1,nspec

  do ip1 = 1,nxgll
  do ip2 = 1,nygll

    xjac2_11 = zero
    xjac2_21 = zero
    xjac2_12 = zero
    xjac2_22 = zero

 do in = 1,ngnod

    nnum = knods(in,ispel)

  xjac2_11 = xjac2_11 + dershape(1,in,ip1,ip2)*coorg(1,nnum)
  xjac2_21 = xjac2_21 + dershape(1,in,ip1,ip2)*coorg(2,nnum)
  xjac2_12 = xjac2_12 + dershape(2,in,ip1,ip2)*coorg(1,nnum)
  xjac2_22 = xjac2_22 + dershape(2,in,ip1,ip2)*coorg(2,nnum)

 enddo

!
!----    invert the jacobian matrix
!

 dvolu(ispel,ip1,ip2) = xjac2_11*xjac2_22 - xjac2_12*xjac2_21

  if (dvolu(ispel,ip1,ip2)  <=  zero) stop 'Error : Jacobian undefined !!'

 xjaci(ispel,1,1,ip1,ip2) =   xjac2_22 / dvolu(ispel,ip1,ip2)
 xjaci(ispel,2,1,ip1,ip2) = - xjac2_21 / dvolu(ispel,ip1,ip2)
 xjaci(ispel,1,2,ip1,ip2) = - xjac2_12 / dvolu(ispel,ip1,ip2)
 xjaci(ispel,2,2,ip1,ip2) =   xjac2_11 / dvolu(ispel,ip1,ip2)

  enddo
  enddo

  enddo

!---- calcul des coordonnees interpolees avec les fonctions de forme
!---- interpolation sur grille reguliere en (xi,eta)

  do i=1,iptsdisp
   xirec(i)  = 2.d0*dble(i-1)/dble(iptsdisp-1) - 1.d0
   etarec(i) = 2.d0*dble(i-1)/dble(iptsdisp-1) - 1.d0
  enddo

!---- calcul des interpolateurs de Lagrange (suppose nxgll = nygll)
  do i=0,nxgll-1
   do k=1,iptsdisp
       flagrange(i,k) = hgll(i,xirec(k),xi,nxgll)
   enddo
  enddo

!
!---- set up the shape functions for the interpolated grid
!
  if(ngnod   ==   4) then
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

  else if(ngnod   ==   9) then
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

  return
  end subroutine q49spec
