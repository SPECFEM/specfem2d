
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.1
!                   ------------------------------
!
!                         Dimitri Komatitsch
!          Universite de Pau et des Pays de l'Adour, France
!
!                          (c) January 2005
!
!========================================================================

  subroutine checkgrid(deltat,f0,t0,initialfield,rsizemin,rsizemax, &
    cpoverdxmin,cpoverdxmax,rlamdaSmin,rlamdaSmax,rlamdaPmin,rlamdaPmax)

!
!----  verification taille des mailles, stabilite et nb de points par lambda
!

  implicit none

  include "constants.h"

  double precision f0,t0
  double precision deltat,rsizemin,rsizemax,cpoverdxmin,cpoverdxmax, &
    rlamdaSmin,rlamdaSmax,rlamdaPmin,rlamdaPmax

  logical initialfield

!
!----  verification taille de grille min et max
!

  print *
  print *,'******************************************'
  print *,'*** Verification parametres simulation ***'
  print *,'******************************************'
  print *
  print *,'*** Taille max grille = ',rsizemax
  print *,'*** Taille min grille = ',rsizemin
  print *,'*** Rapport max/min = ',rsizemax/rsizemin
  print *
  print *,'*** Stabilite max vitesse P = ',cpoverdxmax*deltat
  print *,'*** Stabilite min vitesse P = ',cpoverdxmin*deltat
  print *

  if(.not. initialfield) then

    print *,' Onset time = ',t0
    print *,' Fundamental period = ',1.d0/f0
    print *,' Fundamental frequency = ',f0
    if(t0 <= 1.d0/f0) then
      stop 'Onset time too small'
    else
      print *,' --> onset time ok'
    endif
    print *,'----'
    print *,' Nb pts / lambda P max f0 = ',NGLLX*rlamdaPmax/f0
    print *,' Nb pts / lambda P min f0 = ',NGLLX*rlamdaPmin/f0
    print *,' Nb pts / lambda P max fmax = ',NGLLX*rlamdaPmax/(2.5d0*f0)
    print *,' Nb pts / lambda P min fmax = ',NGLLX*rlamdaPmin/(2.5d0*f0)
    print *,'----'
    print *,' Nb pts / lambda S max f0 = ',NGLLX*rlamdaSmax/f0
    print *,' Nb pts / lambda S min f0 = ',NGLLX*rlamdaSmin/f0
    print *,' Nb pts / lambda S max fmax = ',NGLLX*rlamdaSmax/(2.5d0*f0)
    print *,' Nb pts / lambda S min fmax = ',NGLLX*rlamdaSmin/(2.5d0*f0)
    print *,'----'

  endif

  end subroutine checkgrid

