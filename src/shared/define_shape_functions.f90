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

  subroutine define_shape_functions(shape2D,dershape2D,xi,gamma,ngnod)

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
!                           Local coordinate system : s,t
!
!=======================================================================

  use constants, only: QUARTER,HALF,ONE,TWO,TINYVAL,NDIM

  implicit none

  integer,intent(in) :: ngnod

  double precision,intent(out) :: shape2D(ngnod)
  double precision,intent(out) :: dershape2D(NDIM,ngnod)
  double precision,intent(in) :: xi,gamma

  ! local parameters
  double precision :: s,t,sp,sm,tp,tm,s2,t2,ss,tt,st

!
!---- set up the shape functions and their local derivatives
!
  s  = xi
  t  = gamma

!----    4-node element
  if (ngnod == 4) then
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
  else if (ngnod == 9) then

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
     call stop_the_code('Error: wrong number of control nodes')
  endif

!--- check the shape functions and their derivatives
  ! sum of shape functions should be one
  if (abs(sum(shape2D)-ONE) > TINYVAL) call stop_the_code('Error shape functions')

  ! sum of derivatives of shape functions should be zero
  if (abs(sum(dershape2D(1,:))) > TINYVAL) call stop_the_code('Error deriv xi shape functions')
  if (abs(sum(dershape2D(2,:))) > TINYVAL) call stop_the_code('Error deriv gamma shape functions')

  end subroutine define_shape_functions

