
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
    cpoverdxmax,rlamdaSmin,rlamdaSmax,rlamdaPmin,rlamdaPmax)

!
!----  verification taille des mailles, stabilite et nb de points par lambda
!

  implicit none

  include "constants.h"

  double precision f0,t0
  double precision deltat,rsizemin,rsizemax,cpoverdxmax, &
    rlamdaSmin,rlamdaSmax,rlamdaPmin,rlamdaPmax

  logical initialfield

!
!----  verification taille de grille min et max
!

  write(IOUT,*)
  write(IOUT,*) '******************************************'
  write(IOUT,*) '*** Verification parametres simulation ***'
  write(IOUT,*) '******************************************'
  write(IOUT,*)
  write(IOUT,*) '*** Taille max grille = ',rsizemax
  write(IOUT,*) '*** Taille min grille = ',rsizemin
  write(IOUT,*) '*** Rapport max/min = ',rsizemax/rsizemin
  write(IOUT,*)
  write(IOUT,*) '*** Stabilite max vitesse P = ',cpoverdxmax*deltat
  write(IOUT,*)

  if(.not. initialfield) then

    write(IOUT,*) ' Onset time = ',t0
    write(IOUT,*) ' Fundamental period = ',1.d0/f0
    write(IOUT,*) ' Fundamental frequency = ',f0
    if(t0 <= 1.d0/f0) then
      stop 'Onset time too small'
    else
      write(IOUT,*) ' --> onset time ok'
    endif
    write(IOUT,*) '----'
    write(IOUT,*) ' Nb pts / lambda P max f0 = ',NGLLX*rlamdaPmax/f0
    write(IOUT,*) ' Nb pts / lambda P min f0 = ',NGLLX*rlamdaPmin/f0
    write(IOUT,*) ' Nb pts / lambda P max fmax = ',NGLLX*rlamdaPmax/(2.5d0*f0)
    write(IOUT,*) ' Nb pts / lambda P min fmax = ',NGLLX*rlamdaPmin/(2.5d0*f0)
    write(IOUT,*) '----'
    write(IOUT,*) ' Nb pts / lambda S max f0 = ',NGLLX*rlamdaSmax/f0
    write(IOUT,*) ' Nb pts / lambda S min f0 = ',NGLLX*rlamdaSmin/f0
    write(IOUT,*) ' Nb pts / lambda S max fmax = ',NGLLX*rlamdaSmax/(2.5d0*f0)
    write(IOUT,*) ' Nb pts / lambda S min fmax = ',NGLLX*rlamdaSmin/(2.5d0*f0)
    write(IOUT,*) '----'

  endif

  end subroutine checkgrid

