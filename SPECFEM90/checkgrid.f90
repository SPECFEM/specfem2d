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

  subroutine checkgrid(deltat,gltfu,nltfl,initialfield)

!
!----  verification taille des mailles, stabilite et nb de points par lambda
!

  use verifs
  use spela202

  implicit none

  integer nltfl
  double precision gltfu(20,nltfl)
  double precision deltat
  logical initialfield

  integer n,isource
  double precision f0,t0

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

!
!----  boucle sur toutes les sources
!

  if(.not. initialfield) then

  do n=1,nltfl

!
!----  determiner type de source
!
  isource = nint(gltfu(1,n))
  f0 = gltfu(5,n)
  t0 = gltfu(6,n)

!
!----  utiliser type de source en temps
!
  if(isource == 6) then
      print *,' Source ',n,': Ricker'
      print *,' Onset time = ',t0
      print *,' Fundamental period = ',1.d0/f0
      print *,' Fundamental frequency = ',f0
      if(t0 <= 1.d0/f0) then
        stop 'Onset time too small'
      else
        print *,' --> onset time ok'
      endif
      print *,'----'
      print *,' Nb pts / lambda P max f0 = ',nxgll*rlamdaPmax/f0
      print *,' Nb pts / lambda P min f0 = ',nxgll*rlamdaPmin/f0
      print *,' Nb pts / lambda P max fmax = ',nxgll*rlamdaPmax/(2.5d0*f0)
      print *,' Nb pts / lambda P min fmax = ',nxgll*rlamdaPmin/(2.5d0*f0)
      print *,'----'
      print *,' Nb pts / lambda S max f0 = ',nxgll*rlamdaSmax/f0
      print *,' Nb pts / lambda S min f0 = ',nxgll*rlamdaSmin/f0
      print *,' Nb pts / lambda S max fmax = ',nxgll*rlamdaSmax/(2.5d0*f0)
      print *,' Nb pts / lambda S min fmax = ',nxgll*rlamdaSmin/(2.5d0*f0)
      print *,'----'
  else if(isource == 7) then
      print *,' Source ',n,': dirac **** not checked ****'
  else
      stop 'Unknown type of source'
  endif

  print *

  enddo

  endif

  return
  end subroutine checkgrid
