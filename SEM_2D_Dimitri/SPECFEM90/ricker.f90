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

  double precision function ricker(t,n,gltfu,nltfl)

! calcul du terme temporel de la source pour un Ricker

  use defpi

  implicit none

  integer nltfl,n
  double precision t
  double precision gltfu(20,nltfl)

  double precision f0,t0,factor,a

! parametres pour la source
  f0 = gltfu(5,n)
  t0 = gltfu(6,n)
  factor = gltfu(7,n)

! Ricker
  a = pi*pi*f0*f0
  ricker = - factor * (1.d0-2.d0*a*(t-t0)**2)*exp(-a*(t-t0)**2)

  return
  end function ricker
