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

  subroutine checksource(gltfu,nltfl,deltat,ncycl)

  use verifs

  implicit none

  integer nltfl,ncycl
  double precision deltat

  double precision gltfu(20,nltfl)

  double precision, external :: ricker,dirac
  integer it,n,isource,i,ncycl2,iseuil
  integer icf(1)
  double precision absfreq,cf,cmaxf

! pour spectre de la source (en simple precision pour routine Netlib)
  real, dimension(:), allocatable :: so,ra,rb,wsave
  real azero,valmax

  print *,'Creating gnuplot file for source time functions'

! arrondir ncycl au nombre pair inferieur
  ncycl2 = ncycl
  if(mod(ncycl2,2) /= 0) ncycl2 = ncycl2 - 1

  allocate(so(ncycl2))
  allocate(ra(ncycl2/2))
  allocate(rb(ncycl2/2))
  allocate(wsave(3*ncycl2+15))

  open(unit=11,file='sources',status='unknown')

! boucle sur tous les pas de temps
  do it=1,ncycl2

! boucle sur toutes les sources
  do n=1,nltfl

! determiner type de source
  isource = nint(gltfu(1,n))

! utiliser type de source en temps
  if(isource == 6) then
      gltfu(19,n) = ricker(it*deltat,n,gltfu,nltfl)
  else if(isource == 7) then
      gltfu(19,n) = dirac(it*deltat,n,gltfu,nltfl)
  else
      gltfu(19,n) = 0.d0
  endif

  enddo

  write(11,*) real(it*deltat),(real(gltfu(19,i)),i=1,nltfl)

  enddo

  close(11)

!
! check central frequency by computing the Fourier transform of the source
!

!! DK DK this part suppressed since does not work with range checking
  goto 333

  azero = 0
  n = 1

  do it=1,ncycl2
    so(it)=sngl(ricker(it*deltat,n,gltfu,nltfl))
  enddo

! initialisation pour routine de FFT de Netlib
  call ezffti(ncycl2,wsave)

! appel routine de FFT de Netlib
  call ezfftf(ncycl2,so,azero,ra,rb,wsave)

! prendre le module de l'amplitude spectrale
  ra(:) = sqrt(ra(:)**2 + rb(:)**2)

! determiner la frequence centrale de la source
  icf = maxloc(ra(1:ncycl2/2 - 1))
  cf = icf(1)/(ncycl2*deltat)

! normaliser le spectre d'amplitude
  valmax = ra(icf(1))
  ra(:) = ra(:) / valmax

! determiner la frequence maximale de la source
  iseuil = ncycl2/2 - 1
  do it=icf(1)+1,ncycl2/2 - 1
    if(ra(it) < sngl(valseuil)) then
      iseuil = it
      exit
    endif
  enddo
  cmaxf = iseuil/(ncycl2*deltat)

  print *,'Estimated central freq of the source is ',cf
  print *,'Estimated max freq of the source is ',cmaxf
  print *,'Nyquist frequency for the sampled time function is ',1.d0/(2.d0*deltat)

! sauvegarde du spectre d'amplitude de la source en Hz au format Gnuplot
  open(unit=10,file='spectrum',status='unknown')
  do it=1,ncycl2/2 - 1
    absfreq = it/(ncycl2*deltat)
    if (absfreq <= sngl(freqmaxrep)) write(10,*) sngl(absfreq),ra(it)
  enddo
  close(10)

 333  continue

  deallocate(so)
  deallocate(ra)
  deallocate(rb)
  deallocate(wsave)

  return
  end subroutine checksource
