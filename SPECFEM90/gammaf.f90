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

  double precision function gammaf (x)
!
!=======================================================================
!
!     G a m m a f :
!     -----------
!
!=======================================================================
!
  use defpi

  implicit none

  double precision x

  double precision, parameter :: half=0.5d0,one=1.d0,two=2.d0

  gammaf = one

  if (x == -half) gammaf = -two*dsqrt(pi)
  if (x ==  half) gammaf =  dsqrt(pi)
  if (x ==  one ) gammaf =  one
  if (x ==  two ) gammaf =  one
  if (x ==  1.5d0) gammaf =  dsqrt(pi)/2.d0
  if (x ==  2.5d0) gammaf =  1.5d0*dsqrt(pi)/2.d0
  if (x ==  3.5d0) gammaf =  2.5d0*1.5d0*dsqrt(pi)/2.d0
  if (x ==  3.d0 ) gammaf =  2.d0
  if (x ==  4.d0 ) gammaf = 6.d0
  if (x ==  5.d0 ) gammaf = 24.d0
  if (x ==  6.d0 ) gammaf = 120.d0

  return
  end function gammaf
