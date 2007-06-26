
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

  subroutine define_shape_functions(shape2D,dershape2D,xi,gamma,ngnod)

!=======================================================================
!
!  Set up the shape functions for the subparametric transformation.
!  The routine can handle 4 or 9 control nodes defined as follows:
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

  integer ngnod

  double precision shape2D(ngnod)
  double precision dershape2D(NDIM,ngnod)
  double precision xi,gamma

  double precision s,t,sp,sm,tp,tm,s2,t2,ss,tt,st

!
!---- set up the shape functions and their local derivatives
!
  s  = xi
  t  = gamma

!----    4-node element
  if(ngnod == 4) then
       sp = s + ONE
       sm = s - ONE
       tp = t + ONE
       tm = t - ONE

!----  corner nodes
       shape2D(1) = QUART * sm * tm
       shape2D(2) = - QUART * sp * tm
       shape2D(3) = QUART * sp * tp
       shape2D(4) = - QUART * sm * tp

       dershape2D(1,1) = QUART * tm
       dershape2D(1,2) = - QUART * tm
       dershape2D(1,3) =  QUART * tp
       dershape2D(1,4) = - QUART * tp

       dershape2D(2,1) = QUART * sm
       dershape2D(2,2) = - QUART * sp
       dershape2D(2,3) =  QUART * sp
       dershape2D(2,4) = - QUART * sm

!----    9-node element
  else if(ngnod == 9) then

       sp = s + ONE
       sm = s - ONE
       tp = t + ONE
       tm = t - ONE
       s2 = s * TWO
       t2 = t * TWO
       ss = s * s
       tt = t * t
       st = s * t

!----  corner nodes
       shape2D(1) = QUART * sm * st * tm
       shape2D(2) = QUART * sp * st * tm
       shape2D(3) = QUART * sp * st * tp
       shape2D(4) = QUART * sm * st * tp

       dershape2D(1,1) = QUART * tm * t * (s2 - ONE)
       dershape2D(1,2) = QUART * tm * t * (s2 + ONE)
       dershape2D(1,3) = QUART * tp * t * (s2 + ONE)
       dershape2D(1,4) = QUART * tp * t * (s2 - ONE)

       dershape2D(2,1) = QUART * sm * s * (t2 - ONE)
       dershape2D(2,2) = QUART * sp * s * (t2 - ONE)
       dershape2D(2,3) = QUART * sp * s * (t2 + ONE)
       dershape2D(2,4) = QUART * sm * s * (t2 + ONE)

!----  midside nodes
       shape2D(5) = HALF * tm * t * (ONE - ss)
       shape2D(6) = HALF * sp * s * (ONE - tt)
       shape2D(7) = HALF * tp * t * (ONE - ss)
       shape2D(8) = HALF * sm * s * (ONE - tt)

       dershape2D(1,5) = -ONE  * st * tm
       dershape2D(1,6) =  HALF * (ONE - tt) * (s2 + ONE)
       dershape2D(1,7) = -ONE  * st * tp
       dershape2D(1,8) =  HALF * (ONE - tt) * (s2 - ONE)

       dershape2D(2,5) =  HALF * (ONE - ss) * (t2 - ONE)
       dershape2D(2,6) = -ONE  * st * sp
       dershape2D(2,7) =  HALF * (ONE - ss) * (t2 + ONE)
       dershape2D(2,8) = -ONE  * st * sm

!----  center node
       shape2D(9) = (ONE - ss) * (ONE - tt)

       dershape2D(1,9) = -ONE * s2 * (ONE - tt)
       dershape2D(2,9) = -ONE * t2 * (ONE - ss)

  else
     call exit_MPI('Error: wrong number of control nodes')
  endif

!--- check the shape functions and their derivatives
! sum of shape functions should be one
! sum of derivaticves of shape functions should be zero
  if(abs(sum(shape2D)-ONE) > TINYVAL) call exit_MPI('error shape functions')
  if(abs(sum(dershape2D(1,:))) > TINYVAL) call exit_MPI('error deriv xi shape functions')
  if(abs(sum(dershape2D(2,:))) > TINYVAL) call exit_MPI('error deriv gamma shape functions')

  end subroutine define_shape_functions

