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

  double precision function dirac(t,n,gltfu,nltfl)

! calcul du terme temporel de la source pour un Dirac

  use timeparams

  implicit none

  integer nltfl,n
  double precision t
  double precision gltfu(20,nltfl)

! "largeur" du dirac (fonction triangle) en nb de pas de temps
  integer, parameter :: ilength=4

  double precision t0,factor

! parametres pour la source
  t0 = gltfu(6,n)
  factor = gltfu(7,n)

! Dirac
  if(dabs(t-t0) <= deltat*dble(ilength)/2.d0) then
    if(t <= t0) then
  dirac = - 2.d0*factor*t/(dble(ilength)*deltat) &
    + factor*(2.d0*t0/(dble(ilength)*deltat) - 1.d0)
    else
  dirac = - 2.d0*factor*t/(dble(-ilength)*deltat) &
    + factor*(2.d0*t0/(dble(-ilength)*deltat) - 1.d0)
    endif
  else
  dirac = 0.d0
  endif

  return
  end function dirac
