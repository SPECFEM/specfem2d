
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

  subroutine q49shape(shape,dershape,xi,yi,ngnod,nxgll,nygll,ndime)
!
!=======================================================================
!
!     "q 4 9 s h a p e" : set up the shape functions and their derivatives
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

  integer ngnod,nxgll,nygll,ndime

  double precision shape(ngnod,nxgll,nxgll)
  double precision dershape(ndime,ngnod,nxgll,nxgll)
  double precision xi(nxgll),yi(nygll)

  double precision, parameter :: &
       zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,quart=0.25d0

  integer l1,l2
  double precision s,sp,sm,t,tp,tm,s2,t2,ss,tt,st

  double precision, external :: hgll

!
!-----------------------------------------------------------------------
!

!
!---- set up the shape functions and their local derivatives
!
  if(ngnod   ==   4) then
!
!----    4-noded rectangular element
!
 do l2 = 1,nygll

    t  = yi(l2)

    do l1 = 1,nxgll

       s  = xi(l1)

       sp = s + one
       sm = s - one
       tp = t + one
       tm = t - one

!
!----          corner nodes
!
       shape(1,l1,l2) = quart * sm * tm
       shape(2,l1,l2) = - quart * sp * tm
       shape(3,l1,l2) = quart * sp * tp
       shape(4,l1,l2) = - quart * sm * tp

       dershape(1,1,l1,l2) = quart * tm
       dershape(1,2,l1,l2) = - quart * tm
       dershape(1,3,l1,l2) =  quart * tp
       dershape(1,4,l1,l2) = - quart * tp

       dershape(2,1,l1,l2) = quart * sm
       dershape(2,2,l1,l2) = - quart * sp
       dershape(2,3,l1,l2) =  quart * sp
       dershape(2,4,l1,l2) = - quart * sm

    enddo
 enddo

  else if(ngnod   ==   9) then
!
!----    9-noded rectangular element
!
 do l2 = 1,nygll

    t  = yi(l2)

    do l1 = 1,nxgll

       s  = xi(l1)

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
       shape(1,l1,l2) = quart * sm * st * tm
       shape(2,l1,l2) = quart * sp * st * tm
       shape(3,l1,l2) = quart * sp * st * tp
       shape(4,l1,l2) = quart * sm * st * tp

       dershape(1,1,l1,l2) = quart * tm * t * (s2 - one)
       dershape(1,2,l1,l2) = quart * tm * t * (s2 + one)
       dershape(1,3,l1,l2) = quart * tp * t * (s2 + one)
       dershape(1,4,l1,l2) = quart * tp * t * (s2 - one)

       dershape(2,1,l1,l2) = quart * sm * s * (t2 - one)
       dershape(2,2,l1,l2) = quart * sp * s * (t2 - one)
       dershape(2,3,l1,l2) = quart * sp * s * (t2 + one)
       dershape(2,4,l1,l2) = quart * sm * s * (t2 + one)
!
!----          midside nodes
!
       shape(5,l1,l2) = half * tm * t * (one - ss)
       shape(6,l1,l2) = half * sp * s * (one - tt)
       shape(7,l1,l2) = half * tp * t * (one - ss)
       shape(8,l1,l2) = half * sm * s * (one - tt)

       dershape(1,5,l1,l2) = -one  * st * tm
       dershape(1,6,l1,l2) =  half * (one - tt) * (s2 + one)
       dershape(1,7,l1,l2) = -one  * st * tp
       dershape(1,8,l1,l2) =  half * (one - tt) * (s2 - one)

       dershape(2,5,l1,l2) =  half * (one - ss) * (t2 - one)
       dershape(2,6,l1,l2) = -one  * st * sp
       dershape(2,7,l1,l2) =  half * (one - ss) * (t2 + one)
       dershape(2,8,l1,l2) = -one  * st * sm
!
!----          center node
!
       shape(9,l1,l2) = (one - ss) * (one - tt)

       dershape(1,9,l1,l2) = -one * s2 * (one - tt)
       dershape(2,9,l1,l2) = -one * t2 * (one - ss)

    enddo
 enddo

  else
        stop 'Error : wrong number of control nodes !!'
  endif

  return
  end subroutine q49shape
