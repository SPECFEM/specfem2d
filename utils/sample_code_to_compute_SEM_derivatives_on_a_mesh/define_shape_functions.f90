
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

  subroutine define_shape_functions(shape2D,dershape2D,xi,gamma)

!=======================================================================
!
!  Set up the shape functions for the superparametric transformation
! (i.e. the geometry is defined with lower-order functions than the unknowns
!  of the problem; see for instance Chapter 16 of the finite-element course of
!  Carlos Felippa in Colorado for a discussion about this).
!  The routine can handle 4 or 9 control nodes defined as follows:
!
!                               4 . . . . 7 . . . . 3
!                               .                   .
!                               .         gamma     .
!                               .                   .
!                               8         9  xi     6
!                               .                   .
!                               .                   .
!                               .                   .
!                               1 . . . . 5 . . . . 2
!
!                           Local coordinate system : (xi,gamma)
!
!=======================================================================

  implicit none

  include "constants.h"

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
       shape2D(1) = QUARTER * sm * tm
       shape2D(2) = - QUARTER * sp * tm
       shape2D(3) = QUARTER * sp * tp
       shape2D(4) = - QUARTER * sm * tp

       dershape2D(1,1) = QUARTER * tm
       dershape2D(1,2) = - QUARTER * tm
       dershape2D(1,3) =  QUARTER * tp
       dershape2D(1,4) = - QUARTER * tp

       dershape2D(2,1) = QUARTER * sm
       dershape2D(2,2) = - QUARTER * sp
       dershape2D(2,3) =  QUARTER * sp
       dershape2D(2,4) = - QUARTER * sm

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
       shape2D(1) = QUARTER * sm * st * tm
       shape2D(2) = QUARTER * sp * st * tm
       shape2D(3) = QUARTER * sp * st * tp
       shape2D(4) = QUARTER * sm * st * tp

       dershape2D(1,1) = QUARTER * tm * t * (s2 - ONE)
       dershape2D(1,2) = QUARTER * tm * t * (s2 + ONE)
       dershape2D(1,3) = QUARTER * tp * t * (s2 + ONE)
       dershape2D(1,4) = QUARTER * tp * t * (s2 - ONE)

       dershape2D(2,1) = QUARTER * sm * s * (t2 - ONE)
       dershape2D(2,2) = QUARTER * sp * s * (t2 - ONE)
       dershape2D(2,3) = QUARTER * sp * s * (t2 + ONE)
       dershape2D(2,4) = QUARTER * sm * s * (t2 + ONE)

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
     stop 'Error: wrong number of control nodes'
  endif

!--- check the shape functions and their derivatives
! sum of shape functions should be one
! sum of derivatives of shape functions should be zero
  if(abs(sum(shape2D)-ONE) > TINYVAL) stop 'error shape functions'
  if(abs(sum(dershape2D(1,:))) > TINYVAL) stop 'error deriv xi shape functions'
  if(abs(sum(dershape2D(2,:))) > TINYVAL) stop 'error deriv gamma shape functions'

  end subroutine define_shape_functions

